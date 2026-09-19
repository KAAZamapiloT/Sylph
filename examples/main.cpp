#include<sylph/sh_coefficients.hpp>

#include<iostream>

int main(){
    sylph::SHCoefficients obj(10);

    int a=obj.index(4,-2);
    std::cout<<a<<"\n";
}