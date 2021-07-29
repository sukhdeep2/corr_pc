#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>


int apply(int (*fun)(int,int), int a, int b) {
  return (*fun)(a,b);
}

int add(int a, int b) {return a + b;}
int multiply(int a, int b) {return a * b;}

int main(int argc, const char* argv[]) {
  int added = apply(add, 2, 4);
  int multiplied = apply(multiply, 2, 4);

  std::cout << "added result: " << added << std::endl;
  std::cout << "multiplied result: " << multiplied << std::endl;
}
