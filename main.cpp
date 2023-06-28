#include <iostream>
#include <fstream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <cstring>
#include <lame.h>
#include <complex>
#include <valarray>
using namespace std;
#pragma pack(2)

struct WAVHEADER
{
    char                chunkId[4];                 // Информация о формате файла (RIFF), Содержит символы “RIFF” в ASCII кодировке;
    unsigned long       ChunkSize;                  // Размер без  chunkId[4];
    char                format[4];                  // Формат потоковых данных (WAVE);
    char                subchunk1Id[4];             // Описание параметров WAV-файла (fmt-chunk);
    unsigned long       subchunk1Size;              // Размер подструктуры  subchunk1  (16 байт);
    unsigned short      audioFormat;                  // Аудио формат (PCM = 1);
    unsigned short      nChannels;                  // Количество каналов (Моно = 1, Стерео = 2);
    unsigned long       SamplesRate;              // Частота дискретизации в Гц;
    unsigned long       ByteRate;                   // Кол-во передаваемых байт в секунду воспроизведения;
    unsigned short      blockAlign;                 // Размер сэмпла в байтах 16 бит = 2 байта моно, 32 бита = 4 байта стерео (включая все каналы);
    unsigned short      BitsPerSample;              // Количество бит в сэмпле. Так называемая “глубина” или точность звучания. 8 бит, 16 бит и т.д. /// битов на отсчет
    char                Subchunk2ID[4];             // Символы "Data", начало чанка данных;
    unsigned long     Subchunk2Size;              // Размер области данных в байтах;

};
const double pi=3.14;
const int WINDOW_X = 1080;
const int WINDOW_Y = 720;
const int WINDOW_FPS = 30;
const int TEXTURE_X = 0.8*WINDOW_X;
const int TEXTURE_Y = 0.8*WINDOW_Y;
using Complex=complex<double>;
void FFT( Complex f[], Complex ftilde[], int log2N )                 // Fast Fourier Transform
{
    int N = 1 << log2N;
    // Reorder
    for ( int i = 0; i < N; i++ )
    {
        int ii = 0, x = i;
        for (int j = 0; j < log2N; j++)
        {
            ii <<= 1;
            ii |= ( x & 1 );
            x >>= 1;
        }
        ftilde[ii] = f[i];
    }
    for ( int s = 1; s <= log2N; s++ )
    {
        int m = 1 << s;
        int m2 = m / 2;
        Complex w = 1.0;
        Complex wm = polar( 1.0, -pi / m2 );
        for ( int j = 0; j < m2; j++ )
        {
            for ( int k = j; k < N; k += m )
            {
                Complex t = w * ftilde[k+m2];
                Complex u =     ftilde[k   ];
                ftilde[k   ] = u + t;
                ftilde[k+m2] = u - t;
            }
            w *= wm;
        }
    }
}
void iFFT( Complex ftilde[], Complex f[], int log2N )                // Inverse Fast Fourier Transform
{
    int N = 1 << log2N;
    for ( int m = 0; m < N; m++ ) ftilde[m] = conj( ftilde[m] );      // Apply conjugate (reversed below)
    FFT( ftilde, f, log2N );
    double factor = 1.0 / N;
    for ( int m = 0; m < N; m++ ) f[m] = conj( f[m] ) * factor;
    for ( int m = 0; m < N; m++ ) ftilde[m] = conj( ftilde[m] );      // Only necessary to reinstate ftilde
}

struct WavData{
public:
    int16_t* data;
    long size;

    WavData(){
        data=NULL;
        size=0;
    }
};

void loadWavFile(const char* fname,WavData *ret,WAVHEADER *header,    FILE* wavEQ){
    FILE* file=fopen(fname,"rb");

    if(file==NULL){
        cout<<"error"<<endl;
        exit(0);
    }else{
        cout<<"File is open: "<<fname<<endl;
        fread(header->chunkId,sizeof(char),4,file);
        fwrite(header->chunkId,sizeof(char),4,wavEQ);
        header->chunkId[4]='\0';
        if(!strcmp(header->chunkId,"RIFF")){
            fread(&header->ChunkSize,sizeof(int16_t),2,file);
            fwrite(&header->ChunkSize,sizeof(int16_t),2,wavEQ);
            fread(header->format,sizeof(char),4,file);
            fwrite(header->format,sizeof(char),4,wavEQ);
            header->format[4]='\0';
            if(!strcmp(header->format,"WAVE")){

                fread(header->subchunk1Id,sizeof(char),4,file);
                fwrite(header->subchunk1Id,sizeof(char),4,wavEQ);

                fread(&header->subchunk1Size,sizeof(int16_t),2,file);
                fwrite(&header->subchunk1Size,sizeof(int16_t),2,wavEQ);

                fread(&header->audioFormat,sizeof(int16_t),1,file);
                fwrite(&header->audioFormat,sizeof(int16_t),1,wavEQ);

                fread(&header->nChannels,sizeof(int16_t),1,file);
                fwrite(&header->nChannels,sizeof(int16_t),1,wavEQ);

                fread(&header->SamplesRate,sizeof(int16_t),2,file);
                fwrite(&header->SamplesRate,sizeof(int16_t),2,wavEQ);

                fread(&header->ByteRate,sizeof(int16_t),2,file);
                fwrite(&header->ByteRate,sizeof(int16_t),2,wavEQ);

                fread(&header->blockAlign,sizeof(int16_t),1,file);
                fwrite(&header->blockAlign,sizeof(int16_t),1,wavEQ);

                fread(&header->BitsPerSample,sizeof(int16_t),1,file);
                fwrite(&header->BitsPerSample,sizeof(int16_t),1,wavEQ);

                fread(header->Subchunk2ID,sizeof(char),4,file);
                fwrite(header->Subchunk2ID,sizeof(char),4,wavEQ);

                fread(&header->Subchunk2Size,sizeof(int16_t),2,file);
                fwrite(&header->Subchunk2Size,sizeof(int16_t),2,wavEQ);

                ret->size=header->Subchunk2Size/sizeof(int16_t);
                ret->data=(int16_t*)malloc(header->Subchunk2Size);
                fread(ret->data,sizeof(int16_t),ret->size,file);

            }
            else{
                cout<<"Error: RIFF File but not a wave file\n";
                exit(0);
            }
        }
        else{
            cout<<"ERROR: not a RIFF file\n";
            exit(0);
        }
    }
    fclose(file);
}

void freeSource(WavData* data){
    free(data->data);
}
void scale(double PixelsPerZone,sf::CircleShape *scale0,sf::CircleShape *scale1,sf::CircleShape *scale2,sf::CircleShape *scale3,sf::CircleShape *scale4,sf::CircleShape *scale5,sf::CircleShape *scale6,sf::CircleShape *scale7,sf::CircleShape *scale8,sf::CircleShape *scale9,sf::CircleShape *scale10){
    (*scale0).setPointCount(64);
    (*scale0).setOrigin(3,3);
    (*scale0).setFillColor(sf::Color::Black);
    (*scale0).setPosition((WINDOW_X-TEXTURE_X)/2+0*PixelsPerZone,WINDOW_Y*0.9);
    (*scale1).setPointCount(64);
    (*scale1).setOrigin(3,3);
    (*scale1).setFillColor(sf::Color::Black);
    (*scale1).setPosition((WINDOW_X-TEXTURE_X)/2+1*PixelsPerZone,WINDOW_Y*0.9);
    (*scale2).setPointCount(64);
    (*scale2).setOrigin(3,3);
    (*scale2).setFillColor(sf::Color::Black);
    (*scale2).setPosition((WINDOW_X-TEXTURE_X)/2+2*PixelsPerZone,WINDOW_Y*0.9);
    (*scale3).setPointCount(64);
    (*scale3).setOrigin(3,3);
    (*scale3).setFillColor(sf::Color::Black);
    (*scale3).setPosition((WINDOW_X-TEXTURE_X)/2+3*PixelsPerZone,WINDOW_Y*0.9);
    (*scale4).setPointCount(64);
    (*scale4).setOrigin(3,3);
    (*scale4).setFillColor(sf::Color::Black);
    (*scale4).setPosition((WINDOW_X-TEXTURE_X)/2+4*PixelsPerZone,WINDOW_Y*0.9);
    (*scale5).setPointCount(64);
    (*scale5).setOrigin(3,3);
    (*scale5).setFillColor(sf::Color::Black);
    (*scale5).setPosition((WINDOW_X-TEXTURE_X)/2+5*PixelsPerZone,WINDOW_Y*0.9);
    (*scale6).setPointCount(64);
    (*scale6).setOrigin(3,3);
    (*scale6).setFillColor(sf::Color::Black);
    (*scale6).setPosition((WINDOW_X-TEXTURE_X)/2+6*PixelsPerZone,WINDOW_Y*0.9);
    (*scale7).setPointCount(64);
    (*scale7).setOrigin(3,3);
    (*scale7).setFillColor(sf::Color::Black);
    (*scale7).setPosition((WINDOW_X-TEXTURE_X)/2+7*PixelsPerZone,WINDOW_Y*0.9);
    (*scale8).setPointCount(64);
    (*scale8).setOrigin(3,3);
    (*scale8).setFillColor(sf::Color::Black);
    (*scale8).setPosition((WINDOW_X-TEXTURE_X)/2+8*PixelsPerZone,WINDOW_Y*0.9);
    (*scale9).setPointCount(64);
    (*scale9).setOrigin(3,3);
    (*scale9).setFillColor(sf::Color::Black);
    (*scale9).setPosition((WINDOW_X-TEXTURE_X)/2+9*PixelsPerZone,WINDOW_Y*0.9);
    (*scale10).setPointCount(64);
    (*scale10).setOrigin(3,3);
    (*scale10).setFillColor(sf::Color::Black);
    (*scale10).setPosition((WINDOW_X-TEXTURE_X)/2+10*PixelsPerZone,WINDOW_Y*0.9);
}
void multipeakScale(int counter0,int delta0,int SampleSize,int songsize,Complex arr[],Complex fft[],Complex multipeak[],WavData *song){
    do {
        for (int i=counter0*delta0;i<SampleSize+delta0*counter0;i+=2){
            arr[(i-delta0*counter0)/2]=(*song).data[i];
        }
        FFT(arr, fft, 15);
        if(counter0==0){
            for (int i = 0; i < SampleSize; i++) {
                multipeak[i]=fft[i];
            }
        }

        for (int i = 0; i < SampleSize; i++) {
            if(multipeak[i].real()<fft[i].real()){
                multipeak[i]=fft[i];
            }
        }
        counter0++;
    }while (((SampleSize+delta0*counter0)<songsize)||(counter0<10));
}

WAVHEADER header;
const int SampleSize=32768;
Complex fft[SampleSize];
double fftAverage[SampleSize];
Complex arr[SampleSize];
Complex inverse[SampleSize];

int main() {
    WavData song;
    ofstream out("data.txt");
    ofstream inversedata("complexdata.txt");
    const char* fname="future.wav";
    FILE* wavEQ=fopen("modified.wav","wb");
    loadWavFile(fname,&song,&header,wavEQ);
    for(long i=0;i<song.size;i+=2){
        out<<song.data[i]<<',';
    }
    out.close();
    sf::SoundBuffer x;
    x.loadFromFile(fname);
    cout<<"SampleCount: "<<x.getSampleCount()<<endl;
    cout<<"SampleRate: "<<x.getSampleRate()<<endl;
    cout<<"Duration: "<<x.getDuration().asSeconds()<<endl;
    int songsize=song.size;
    ofstream fftData("fft.txt");
    sf::Music music;
    music.openFromFile(fname);

    sf::RenderWindow window(sf::VideoMode(1080, 720), "Music Visualiser");
    window.setFramerateLimit(WINDOW_FPS);
    sf::Event event;
    sf::RectangleShape timeline(sf::Vector2f(TEXTURE_X, 1));
    timeline.setFillColor(sf::Color::Black);
    timeline.setPosition((WINDOW_X-TEXTURE_X)/2, 0.9*WINDOW_Y);
    //Seek
    sf::CircleShape seek(3);
    seek.setPointCount(64);
    seek.setOrigin(3,3);
    seek.setFillColor(sf::Color::Red);
    //Graph components
    sf::Vertex vertices[TEXTURE_X];
    for (int i=0; i<TEXTURE_X; i++){
        vertices[i].position.y = TEXTURE_Y/2;
        vertices[i].color = sf::Color(0,0,0,255);
    }
    sf::VertexBuffer vb(sf::LineStrip);
    vb.create(TEXTURE_X);
    sf::RenderTexture rendergraph;
    rendergraph.create(TEXTURE_X, TEXTURE_Y);
    sf::Sprite graph;
    graph.setTexture(rendergraph.getTexture());
    graph.setPosition((WINDOW_X-TEXTURE_X)/2, (WINDOW_Y-TEXTURE_Y)*0.2);

    sf::CircleShape point(2.f);
    point.setFillColor(sf::Color::Red);
    int dur = x.getDuration().asSeconds();
    /*int counter0=0;
    int delta0=header.nChannels*(header.SamplesRate/WINDOW_FPS);
    Complex multipeak[SampleSize];
    multipeakScale(counter0, delta0,SampleSize,songsize,arr,fft,multipeak,&song);
    for (int i=0;i<SampleSize;i++){
        cout<<multipeak[i].real()<<endl;
    }*/
    int sc[11];
    sc[0]=0;
    sc[1]=20;
    sc[2]=50;
    sc[3]=100;
    sc[4]=200;
    sc[5]=500;
    sc[6]=1000;
    sc[7]=2000;
    sc[8]=5000;
    sc[9]=10000;
    sc[10]=20000;
    double PixelsPerZone=(double)TEXTURE_X/10;
    sf::CircleShape scale0(3),scale1(3),scale2(3),scale3(3),scale4(3),scale5(3),scale6(3),scale7(3),scale8(3),scale9(3),scale10(3);
    scale(PixelsPerZone,&scale0,&scale1,&scale2,&scale3,&scale4,&scale5,&scale6,&scale7,&scale8,&scale9,&scale10);
    sf::Font font;
    font.loadFromFile("calibri.ttf");
    sf::Text txt[11];
    sf::Text txtSec;
    for (int i = 0; i < 11; i++) {
        txt[i].setFont(font);
        txt[i].setCharacterSize(20);
        txt[i].setFillColor(sf::Color::Black);
        txt[i].setStyle(sf::Text::Regular);
        txt[i].setPosition((WINDOW_X-TEXTURE_X)/2+i*PixelsPerZone,WINDOW_Y*0.85);
    }
    txt[0].setString("0");
    txt[1].setString("20");
    txt[2].setString("50");
    txt[3].setString("100");
    txt[4].setString("200");
    txt[5].setString("500");
    txt[6].setString("1k");
    txt[7].setString("2k");
    txt[8].setString("5k");
    txt[9].setString("10k");
    txt[10].setString("20k");
    txtSec.setFont(font);
    txtSec.setFillColor(sf::Color::Black);
    txtSec.setStyle(sf::Text::Regular);
    txtSec.setPosition((WINDOW_X-TEXTURE_X)/2,WINDOW_Y*0.95);
    double koef1=1,koef2=1,koef3=1,koef4=1;
    sf::Clock clock;
    while(window.isOpen()){
        //music.play();
        while(window.pollEvent(event)) {
            if (event.type == 0 || (event.type == 5 && event.key.code == 36))
                window.close();
            else if (event.type == 5 && event.key.code == 57) {
                if (music.getStatus() == 2)
                    music.pause();
                else {
                    music.play();
                }

            }
        }
        int delta=header.nChannels*(header.SamplesRate/WINDOW_FPS);

        float time1;
        int f=(songsize-SampleSize)/delta;
        int counter=0;
        int counter3=0;
        int flag=0;
        int flag1=0;
        sf::Clock clock2;
        do {
            window.clear(sf::Color::White);
            rendergraph.clear(sf::Color::White);
            while(window.pollEvent(event)) {
                if (event.type == sf::Event::Closed || (event.type == 5 && event.key.code == 36)){
                    window.close();
                    music.stop();
                    flag1=1;
                }else{
                    if (event.type == 5 && event.key.code == 57) {
                        if (music.getStatus() == 2){
                            music.pause();
                        }else {
                            music.play();
                        }
                    }
                }
            }
            if(flag1){
                break;
            }
            for (int i=counter*delta;i<2*SampleSize+delta*counter;i+=2){
                if(i>songsize){
                    flag=1;
                    break;
                }
                arr[(i-delta*counter)/2].real(song.data[i]);
                arr[(i-delta*counter)/2].imag(0.0);
            }

            FFT(arr, fft, 15);
            counter++;
            int n=2;
            for (int j = 0; j < 10*n; j++){
                fftAverage[j]=0;
            }
            for(int i=0; i<TEXTURE_X; i++){
                for (int j = 0; j < 10*n; j++) {
                    if((PixelsPerZone/n*(double)j<=(double )i)&&(PixelsPerZone/n*((double)j+1)>=(double )i)){
                        double step=(sc[(j+1)/n]-sc[j/n])/PixelsPerZone/2*(i-PixelsPerZone*j/n);
                        fftAverage[j]+=2.00/SampleSize*sqrt(pow(fft[(int)(sc[j/n]/2+step)].real(),2)+pow(fft[(int)(sc[j/n]/2+step)].imag(),2))/(PixelsPerZone/n);
                        break;
                    }
                }
            }
            for(int i=0; i<TEXTURE_X; i++){
                for (int j = 0; j < 10*n; j++) {
                    if((PixelsPerZone*j/n<=(double )i)&&(PixelsPerZone*(j+1)/n>=(double )i)){
                        if(j>=5*n){
                            if(j>=8*n){
                                if(fftAverage[j]*koef1>TEXTURE_Y*3/4){
                                    koef1-=0.1*koef1;
                                }
                                vertices[i].position = sf::Vector2f(i, fftAverage[j]*koef1);
                            }else{
                                if(fftAverage[j]*koef2>TEXTURE_Y*3/4){
                                    koef2-=0.1*koef2;
                                }
                                vertices[i].position = sf::Vector2f(i, fftAverage[j]*koef2);
                            }
                        }else{
                            if(j<2*n){
                                if(fftAverage[j]*koef3>TEXTURE_Y*3/4){
                                    koef3-=0.1*koef3;
                                }
                                vertices[i].position = sf::Vector2f(i, fftAverage[j]*koef3);
                            }else{
                                if(fftAverage[j]*koef4>TEXTURE_Y*3/4){
                                    koef4-=0.1*koef4;
                                }
                                vertices[i].position = sf::Vector2f(i, fftAverage[j]*koef4);

                            }
                        }
                        break;
                    }
                }
            }

            vb.update(vertices);
            int now_sec = music.getPlayingOffset().asSeconds();
            int pos = (WINDOW_X-TEXTURE_X)/2+now_sec*TEXTURE_X/dur;
            seek.setPosition(pos,WINDOW_Y*0.9);

            rendergraph.draw(vb);
            window.draw(graph);
            window.draw(timeline);
            window.draw(seek);
            window.draw(scale0);
            window.draw(scale1);
            window.draw(scale2);
            window.draw(scale3);
            window.draw(scale4);
            window.draw(scale5);
            window.draw(scale6);
            window.draw(scale7);
            window.draw(scale8);
            window.draw(scale9);
            window.draw(scale10);
            for (int i = 0; i < 11; i++) {
                window.draw(txt[i]);
            }
            for (int i = 0; i < SampleSize; i++) {
                if((i<500/2)||(i>(SampleSize-500/2))){
                    fft[i].real(0);
                    fft[i].imag(0);
                }
            }
            if((counter-1)==0 ){
                iFFT( fft,inverse, 15);
                for (int i = 0; i < SampleSize/2; i++) {
                    short int v=(short int)inverse[i].real();
                    fwrite(&v, sizeof(short int), 1, wavEQ);
                    fwrite(&v, sizeof(short int), 1, wavEQ);
                }
            }else {
                if (flag != 1) {
                    iFFT(fft, inverse, 15);
                    for (int i = SampleSize / 2 - delta / 2; i < SampleSize / 2; i++) {
                        short int v = (short int) inverse[i].real();
                        fwrite(&v, sizeof(short int), 1, wavEQ);
                        fwrite(&v, sizeof(short int), 1, wavEQ);
                    }
                } else {
                    iFFT(fft, inverse, 15);
                    for (int i = SampleSize / 2 - delta / 2; i < SampleSize; i++) {
                        //short int v=(short int)sqrt(pow(inverse[i].real(),2)+pow(inverse[i].imag(),2));
                        short int v = (short int) inverse[i].real();
                        fwrite(&v, sizeof(short int), 1, wavEQ);
                        fwrite(&v, sizeof(short int), 1, wavEQ);
                    }
                }
            }
            int f=counter/30;
            string textS;
            textS=to_string(f);
            txtSec.setString(textS);
            window.draw(txtSec);
            if(counter==1){
                music.play();
                clock.restart();
            }
            time1=clock.getElapsedTime().asMicroseconds();
            if(time1<33333){
                sf::sleep(sf::microseconds(33333-time1));
            }else{
                if(time1>35000){
                    music.pause();
                    sf::sleep(sf::microseconds(time1-35000));
                    music.play();
                }
            }
            window.display();
            clock.restart();
        }while ((delta*counter)<songsize);
        break;
    }
    fftData.close();
    inversedata.close();
    freeSource(&song);
    return 0;
}
