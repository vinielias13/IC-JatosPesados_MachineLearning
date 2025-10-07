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


bool traceMotherForBottomHadron(TClonesArray *particles, Int_t index, Int_t &idxLastCharmQuark);

Int_t idxlast = -1;

// Declaração global
Double_t phi_last_quark, eta_last_quark, pt_last_quark;

Int_t findLastQuark(TClonesArray *particulas, Int_t index = -1);

void DadosCSVBBBar()
{

    std::ofstream file("NovoBBar_Tagging.csv");
    file << "PtJet,nConstituent,Angularidade,SigmaKT,MaxRho,N_Vertices_S,FracLeadPt,Label" << std::endl;

    TFile *rootfile = TFile::Open("BTagging.root", "READ");

    if (!rootfile || !rootfile->IsOpen())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    TTree *ttree = dynamic_cast<TTree *>(rootfile->Get("tree"));

    if (!ttree)
    {
        std::cerr << "Error: Could not find the TTree in the file." << std::endl;
        rootfile->Close();
        return;
    }


    // O Array com as variaveis de interesse
    TClonesArray *particulas = new TClonesArray("MyPart");
    TClonesArray *jatos = new TClonesArray("MyJet");
    TClonesArray *constituintes = new TClonesArray("MyConstituent");

    // Coluna com as variaveis
    ttree->SetBranchAddress("particulas", &particulas);
    ttree->SetBranchAddress("jatos", &jatos);
    ttree->SetBranchAddress("constituintes", &constituintes);

    Long64_t nev = ttree->GetEntries();

    // Lista de PDG codes para os hádrons de B que serão aceitos
    const int nBottomHadrons = 8;
    int BottomHadrons[nBottomHadrons] = {511, 521, 513, 523, 531, 533, 541, 543};

    //---------------------------------------------------------------------------------------------------------
    // Loop sobre os eventos
    //---------------------------------------------------------------------------------------------------------

    for (Long64_t i = 0; i < 100000; i++)
    {

        ttree->GetEntry(i);
        Int_t start = 0;

        //---------------------------------------------------------------------------------------------------------
        // Loop sobre os jatos
        //---------------------------------------------------------------------------------------------------------

        for (int ijet = 0; ijet < jatos->GetEntries(); ijet++)
        {

            // Inicializacao de variaveis e reinicializacao a cada evento
            Double_t Ang = 0, sigmaKT = 0, somaPtConstituente = 0, Rho = 0, MaxRho = 0;
            Int_t NSecV = 0;
            int candidateIndex = -1;

            // Boleano que se torna verdadeiro se o jato for proveninte de quark bottom
            bool jet_has_bottom = false;

            MyJet *jet = static_cast<MyJet *>(jatos->At(ijet));

            Double_t jetPt = jet->jetPt;
            Double_t jetEta = jet->jetEta;
            Double_t jetPhi = jet->jetPhi;
            Double_t jetMass = jet->jetMass;
            Double_t jetPx = jet->jetPx;
            Double_t jetPy = jet->jetPy;
            Double_t jetPz = jet->jetPz;
            Double_t jetE = jet->jetE;
            Int_t nConstituents = jet->nConstituent;
            Int_t Ipjet = jet->Ipjet;

            Double_t soma = 0;

            Double_t pTLeadConstituent = 0; // pT do constituinte líder

            Int_t idxQuarkAntesHadron = -1;

            Double_t ptQuark = -1;

            for (Int_t iconst = 0; iconst < (nConstituents); iconst++)
            {

                Int_t idx = start + iconst;
                MyConstituent *constituent = static_cast<MyConstituent *>(constituintes->At(idx));

                Double_t pt = constituent->pt;
                Double_t px = constituent->px;
                Double_t py = constituent->py;
                Double_t pz = constituent->pz;
                Double_t eta = constituent->eta;
                Double_t phi = constituent->phi;
                Double_t E = constituent->E;

                // soma += constituent->pt;;

                Double_t vx = constituent->vx;
                Double_t vy = constituent->vy;
                Double_t vz = constituent->vz;

                Int_t pdg = constituent->pdg;
                Int_t motherIndexconst = constituent->motherIndex;
                Int_t ip = constituent->ipTclones;

                if (!jet_has_bottom)
                {

                    //---------------------------------------------------------------------------------------------------------
                    // Loop das maes
                    //---------------------------------------------------------------------------------------------------------

                    // Recupera particulas no TClonesArrey
                    MyPart *mp = static_cast<MyPart *>(particulas->At(ip));

                    // Indice da primeira mae do constituinte
                    Int_t motherIndex = mp->fMotherIndex;


                    while (motherIndex >= 7)
                    {

                        MyPart *mom = static_cast<MyPart *>(particulas->At(motherIndex));

                        // Obten o valor absoluto do PDG
                        int absMotherPdg = TMath::Abs(mom->fPdg);

                        // Se achamos um hádron de bottom na árvore genealógica, seguimos até o quark bottom original
                        for (int i = 0; i < nBottomHadrons; i++)
                        {

                            // Se alguma mae estiver na lista dos hadrons de bottom entra no If
                            if (absMotherPdg == BottomHadrons[i])
                            {
                                MyPart *mp = static_cast<MyPart *>(particulas->At(ip));

                                // pega o pdg diretamente de mp
                                int pdg = mp->fPdg;

                                // Int_t idxQuarkAntesHadron = -1;
                                // Se a funcao retornar True sabeberemos que o jato e de bottom
                                Int_t idxLastCharmQuark = -1;
                                if (traceMotherForBottomHadron(particulas, motherIndex, idxQuarkAntesHadron))
                                {

                                    // Para utilizar na construcao do match
                                    idxlast = idxQuarkAntesHadron;
                                    jet_has_bottom = true;
                                    break; // Se econtrarmos pelo menos um hadron de charm saimos do loop
                                }

                            }

                        }

                        if (jet_has_bottom)
                        {

                            break; // isso sai do while
                        }

                        // Atualiza o index e constinua subindo as cadeias genealogicas
                        motherIndex = mom->fMotherIndex;

                    } // Fim do loop das maes
                }

                //---------------------------------------------------------------------------------------------------------
                // Variaveis derivadas utilizadas no aprendizado de maquinas
                //---------------------------------------------------------------------------------------------------------

                //---------------------------------------------------------------------------------------------------------
                // Cálculo de Rho (distância radial no plano XY)
                Double_t Rho = TMath::Sqrt(pow(vx, 2) + pow(vy, 2));


                if ((Rho >= 0.01) && (Rho <= 1))
                {
                    NSecV++;
                    if ((Rho <= 1) && (Rho > MaxRho)){
                        MaxRho = Rho;
                    }
                }
                //---------------------------------------------------------------------------------------------------------

                //---------------------------------------------------------------------------------------------------------
                // Definindo um Phi entre -pi e pi
                Double_t deltaPhi = TMath::Abs(phi - jetPhi);

                if (deltaPhi > TMath::Pi())
                {
                    deltaPhi = 2 * TMath::Pi() - deltaPhi;
                }

                // Calculo da distancia angular R entre constituinte e jato
                Double_t R = TMath::Sqrt(TMath::Power(eta - jetEta, 2) + TMath::Power(deltaPhi, 2));
                //---------------------------------------------------------------------------------------------------------

                //---------------------------------------------------------------------------------------------------------
                // Se condicao aceita vou atualizar o pT do lead
                if (pTLeadConstituent < pt)
                {
                    pTLeadConstituent = pt;
                }
                //---------------------------------------------------------------------------------------------------------
                
                // Soma da Angulatidade e dispersao Sigma
                Ang += R;
                sigmaKT += TMath::Power(pt - (jetPt / nConstituents), 2);

            

            } // Fim do loop dos constituintes

            start += nConstituents;

            //---------------------------------------------------------------------------------------------------------

            Double_t fracLeadPt = pTLeadConstituent / jetPt;

            Double_t Angularidade = Ang / nConstituents;

            Double_t DesvioSigmaKt = TMath::Sqrt(sigmaKT/nConstituents);


            //---------------------------------------------------------------------------------------------------------

            if ((jetPt > 5))
            {

                if (jet_has_bottom)
                {
                    if (idxlast >= 0)
                    {
                        MyPart *quarkBottom = static_cast<MyPart *>(particulas->At(idxlast));

                        double_t ptQuark = quarkBottom->fPt;
                        Double_t PhiQuark = quarkBottom->fPhi;
                        Double_t EtaQuark = quarkBottom->fEta;

                        if (ptQuark < 0)
                        {
                            std::cerr << "ptQuark inválido!" << std::endl;
                        }

                        // Definindo um Phi entre -pi e pi
                        Double_t deltaPhi = TMath::Abs(PhiQuark - jetPhi);

                        if (deltaPhi > TMath::Pi())
                        {
                            deltaPhi = 2 * TMath::Pi() - deltaPhi;
                        }

                        // Calcular a distância
                        Double_t distancia = TMath::Sqrt(TMath::Power(deltaPhi, 2) + TMath::Power(EtaQuark - jetEta, 2));


                        if ((distancia <= 0.20))
                        {
                            
                            Int_t Label = 1;
                            file << jetPt << "," << nConstituents << "," << Angularidade << "," << DesvioSigmaKt << "," << MaxRho << "," << NSecV << "," << fracLeadPt << "," << Label << std::endl;
                        
                        }
                    }

                        
                }
                else
                {
                    // Caso nao seja um jato de Bottom
                    Int_t Label = 0;
                    file << jetPt << "," << nConstituents << "," << Angularidade << "," << DesvioSigmaKt << "," << MaxRho << "," << NSecV << "," << fracLeadPt << "," << Label << std::endl;

                }
            }

        } // Fim do loop dos jatos

    } // Fim do loop dos eventos

    file.close();

}

bool traceMotherForBottomHadron(TClonesArray *particulas, Int_t idx, Int_t &idxLastCharmQuark)
{
    MyPart *mp = static_cast<MyPart *>(particulas->At(idx));
    if (!mp)
        return false;

    Int_t momIdx = mp->fMotherIndex;
    if (momIdx < 0)
        return false;

    MyPart *mom = static_cast<MyPart *>(particulas->At(momIdx));
    if (!mom)
        return false;

    if ((idxLastCharmQuark == -1) && (TMath::Abs(mom->fPdg) == 5))
    {
        idxLastCharmQuark = momIdx;
    }

    // Parou em um quark Bottom
    if ((momIdx == 4 || momIdx == 5) && (TMath::Abs(mom->fPdg) == 5))
    {

        return true;
    }

    // Senão, sobe recursivamente
    return traceMotherForBottomHadron(particulas, momIdx, idxLastCharmQuark);
}
