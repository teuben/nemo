#include <iostream>
#include <exception.h>

using namespace WDutils;

int main()
{
  std::cout<<" Hello World\n";
  std::cout<<"\"RunInfo::host()\" = "<<RunInfo::host()<<'\n';
  std::cout<<"\"RunInfo::name()\" = "<<RunInfo::name()<<'\n';
  std::cout<<"\"RunInfo::cmd ()\" = "<<RunInfo::cmd ()<<'\n';
}
