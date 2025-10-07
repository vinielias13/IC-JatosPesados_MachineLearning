#ifndef MyConstituent_h
#define MyConstituent_h

#include <TObject.h>

class MyConstituent : public TObject
{
 public:
  MyConstituent() :
    TObject(), pt(0), px(0), py(0), pz(0),
    eta(0), phi(0), E(0),
    vx(0), vy(0), vz(0),
    pdg(0), motherIndex(-1), ipTclones(0) {}

  // Momento e energia
  Double32_t pt;
  Double32_t px;
  Double32_t py;
  Double32_t pz;
  Double32_t eta;
  Double32_t phi;
  Double32_t E;

  // Informação do vértice
  Double32_t vx;
  Double32_t vy;
  Double32_t vz;

  // Informações adicionais
  Int_t pdg;
  Int_t motherIndex;
  Int_t ipTclones;

  ClassDef(MyConstituent, 1)
};

#endif
