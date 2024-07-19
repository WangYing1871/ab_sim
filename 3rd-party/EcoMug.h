/////////////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                         //
// Copyright (C) 2021 Davide Pagano <davide.pagano@unibs.it>                       //
// EcoMug is based on the following work:                                          //
// D. Pagano, G. Bonomi, A. Donzella, A. Zenoni, G. Zumerle, N. Zurlo,             //
// "EcoMug: an Efficient COsmic MUon Generator for cosmic-ray muons applications", //
// doi:10.1016/j.nima.2021.165732                                                  //
//                                                                                 //
// This program is free software: you can redistribute it and/or modify            //
// it under the terms of the GNU General Public License as published by            //
// the Free Software Foundation, either version 3 of the License, or               //
// (at your option) any later version.                                             //
//                                                                                 //
// This program is distributed in the hope that it will be useful,                 //
// but WITHOUT ANY WARRANTY; without even the implied warranty of                  //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   //
// GNU General Public License for more details.                                    //
//                                                                                 //
// You should have received a copy of the GNU General Public License               //
// along with this program.  If not, see <https://www.gnu.org/licenses/>.          //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef EcoMug_H
#define EcoMug_H 1 
#include <cmath>
#include <array>
#include <random>
#include <functional>

//! Fast generation of random numbers
//! This class is based on the xoroshiro128+ generator.
//! https://prng.di.unimi.it/
struct EMRandom{
  EMRandom(){
    std::random_device rd;
    std::mt19937 gen(std::invoke(rd));
    std::uniform_int_distribution<uint64_t> 
      dis(0,std::numeric_limits<uint64_t>::max());
    m_s[0] = dis(gen); m_s[1] = dis(gen); }

  void SetSeed(uint64_t v){ m_s[0] = v; m_s[1] = v; }
  double GenerateRandomDouble(){ return to_double(next()); }
  double GenerateRandomDouble(double x1, double x2){ return (x2-x1)*GenerateRandomDouble()+x1; }
  int64_t rotl(uint64_t const x, int k){ return (x<<k)|(x>>(64-k)); }
  uint64_t next(){
    uint64_t s0 = m_s[0]; uint64_t s1 = m_s[1];
    std::uint64_t const result = s0+s1;
    s1 ^= s0;
    m_s[0] = rotl(s0,55) ^ s1 ^ (s1<<14);
    m_s[1] = rotl(s1,36); }
  double to_double(uint64_t x){
    uint64_t rt = ((std::uint64_t)(0x3FF)<<52)|(x>>12);
    return *reinterpret_cast<double*>(std::addressof(rt))-1.0;}



public:
  uint64_t m_s[2];
};
//---------------------------------------------------------------------
class EMMaximization{
  
};
#endif/! https://prng.di.unimi.it/
