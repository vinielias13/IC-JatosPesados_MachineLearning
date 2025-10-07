#include "iostream"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "MyPart.h"
#include "MyJet.h"
#include "MyConstituent.h"
#include "TFile.h"


using namespace std;


void CreatQCD(Int_t nev = 500000, Int_t ndeb = 1)
{
    //caregando bibiliotecas do pythia
    gSystem->Load("LibEG");
    gSystem->Load("LibEGPythia8");

    //Criando arquivo e TTree 
    TFile *outfile = TFile::Open("HardQCDO.root", "RECREATE");
    
    //Criando TTree 
    TTree *tree = new TTree("tree", "HardQCDO Trees");

    TClonesArray *particulas = new TClonesArray("MyPart");
    TClonesArray *jatos =  new TClonesArray("MyJet");
    TClonesArray *constituintes = new TClonesArray("MyConstituent");

    //Criando galho para particulas
    tree->Branch("particulas", particulas, 0, 0);
    tree->Branch("jatos", jatos, 0, 0);
    tree->Branch("constituintes", constituintes, 0, 0);

    
    //---------------------------------------------------------------------------------------------------------
    // configuracoes e inicializacoes do Pythia8:
    //---------------------------------------------------------------------------------------------------------

    // Cria o TClonesArray
    TClonesArray *particles = new TClonesArray("TParticle", 1000);

    // Cria o TPythia8
    TPythia8 *pythia = new TPythia8();

    pythia->ReadString("HardQCD:all = on");
    pythia->ReadString("Random:setSeed = on");
    pythia->ReadString("Random:seed = 42");
    pythia->Initialize(2212, 2212, 14000);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacoes e configuracoes do FastJet:
    //---------------------------------------------------------------------------------------------------------

    Double_t R = 0.4;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

    std::vector<fastjet::PseudoJet> fjParticles;
   

    //---------------------------------------------------------------------------------------------------------
    // Loop sobre os eventos
    //---------------------------------------------------------------------------------------------------------

    for (Int_t iev = 0; iev < nev; iev++){
        
   // cout << "//---------------------------------------------------------------------------------------------------------//" <<endl;
        
        particulas->Clear();
        jatos->Clear();
        constituintes->Clear();
        fjParticles.clear();


        // Gerando evento
        pythia->GenerateEvent(); 

        // Importando partículas do evento
        pythia->ImportParticles(particles, "All");
        
        if (iev < ndeb){ 

            // Listando eventos para depuração
            pythia->EventListing(); 
        
        }

        Int_t np = particles->GetEntriesFast();
        Int_t nfp = 0;
        Int_t nfp2 = 0;
        Int_t nfp3 = 0;
        

        //---------------------------------------------------------------------------------------------------------
        // Loop sobre as particulas
        //---------------------------------------------------------------------------------------------------------

        for (Int_t ip = 0; ip < np; ip++) 
        {
            
            // Resgatando particulas e Pegando o estatus, e pdg
            TParticle *part = (TParticle*)particles->At(ip);
            //Int_t pdg = part->GetPdgCode();

            Int_t ist = part->GetStatusCode();
            Int_t pdg = part->GetPdgCode();

        

            MyPart *mp = static_cast<MyPart *>(particulas->New(nfp++));

            // if (ist<=0) continue;
            
            mp->fPt = part->Pt();               
            mp->fEta = part->Eta();              
            mp->fPhi = part->Phi();                          
            mp->fPx = part->Px();               
            mp->fPy = part->Py();               
            mp->fPz = part->Pz();               
            mp->fE = part->Energy();                
            mp->fPdg = pdg;              
            mp->fIp = ip;
            mp->fStatus = ist;
            mp->fMotherIndex   = part->GetFirstMother();

           

            if ((ist > 0) && (TMath::Abs(part->Eta()) < 2))  {
                fastjet::PseudoJet pj(part->Px(), part->Py(), part->Pz(), part->Energy());
                pj.set_user_index(ip); // Armazena o índice da partícula original
                fjParticles.push_back(pj);
            }

          
        }

        // Clusteriza os jatos e obtensao da lista
        fastjet::ClusterSequence clust_seq(fjParticles, jet_def);
        std::vector<fastjet::PseudoJet> jets = clust_seq.inclusive_jets();

        //---------------------------------------------------------------------------------------------------------
        // Loop sobre os jatos
        //---------------------------------------------------------------------------------------------------------
        
        for (int ijet=0; ijet<jets.size(); ijet++) 
        {   

            fastjet::PseudoJet jet=jets[ijet];

            MyJet *mj = static_cast<MyJet *>(jatos->New(nfp2++));

            mj->jetPt   = jet.pt();
            mj->jetEta  = jet.eta();
            mj->jetPhi  = jet.phi_std();
            mj->jetMass = jet.m();
            mj->jetPx   = jet.px();
            mj->jetPy   = jet.py();
            mj->jetPz   = jet.pz();
            mj->jetE    = jet.E();
            mj->nConstituent = jet.constituents().size();
            mj->Ipjet = ijet;

            Double_t soma = 0;

            // Obtém os constituintes do jato
            std::vector<fastjet::PseudoJet> constituents = jet.constituents();

            // Obtém os constituintes do jato
            Int_t nConstituent = constituents.size();  // Número de constituintes   


            //---------------------------------------------------------------------------------------------------------
            // Loop sobre os constituintes
            //---------------------------------------------------------------------------------------------------------

            for (Int_t iconst=0; iconst < nConstituent; iconst++) 
            {

                fastjet::PseudoJet constituent = constituents[iconst];

                Int_t ip = constituent.user_index(); // Índice da partícula original
               
                TParticle *part = (TParticle *)particles->At(ip);
                
                MyConstituent *mc = static_cast<MyConstituent *>(constituintes->New(nfp3++));

                mc->pt  = constituent.pt();
                mc->px  = constituent.px();
                mc->py  = constituent.py();
                mc->pz  = constituent.pz();
                mc->eta = constituent.eta();
                mc->phi = constituent.phi();
                mc->E   = constituent.e();

                soma += constituent.pt();
        
                mc->vx = part->Vx();  // posição x do vértice
                mc->vy = part->Vy(); 
                mc->vz = part->Vz();  
                mc->pdg = part->GetPdgCode();
                mc->motherIndex = part->GetFirstMother();
                mc->ipTclones = ip;
            }


        }
  
        tree->Fill();    
    }     

    pythia->PrintStatistics();
   
    tree->Write();
    //outfile->Close();

    // Limpeza da memória
    delete particulas;
    delete jatos;
    delete constituintes;
    delete particles;
    delete pythia;

}