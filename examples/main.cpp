#include<iostream>

#include<gl/GL.h>
#include<SDL3/SDL.h>


int main(){

    int w=500;
    int h=500;

SDL_Window*window=SDL_CreateWindow("SDL3 + MSVC",
        w, h,
        SDL_WINDOW_RESIZABLE);

        
    
    return 0;
}