#ifndef MyPart_h
#define MyPart_h

class MyPart;

class MyPart: public TObject
{
 public: MyPart() :
  TObject(), fPt(0), fEta(0), fPhi(0), fMass(0), fPx(0), fPy(0), fPz(0), fE(0),  fPdg(0), fIp(-1), fStatus(0), eta_last_quark(0), phi_last_quark(0), pt_last_quark(0),
  fMotherIndex(-1), fDaughter1Index(-1), fDaughter2Index(-1) {;}

 public:
  Double32_t    fPt;               //[0,0,16] pt
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fMass;             //[0,0,16] mass
  Double32_t    fPx;               //[0,0,16] px
  Double32_t    fPy;               //[0,0,16] py
  Double32_t    fPz;               //[0,0,16] pz
  Double32_t    fE;                //[0,0,16] energy
  Int_t         fPdg;              //pdg code
  Int_t         fIp;
  Int_t         fStatus;

  Double32_t    eta_last_quark;
  Double32_t    phi_last_quark;
  Double32_t    pt_last_quark;

  Int_t         fMotherIndex;      // Índice da partícula mãe
  Int_t         fDaughter1Index;   // Índice da primeira filha
  Int_t         fDaughter2Index;   // Índice da segunda filha
  
  ClassDef(MyPart, 1)
    };

#endif