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

// variaveis globais do hadron aceito
Double_t phi_charmHadron, eta_charmHadron, pt_charmHadron;

bool traceMotherForBottomHadron(TClonesArray *particles, Int_t index, Int_t &idxLastCharmQuark);

Int_t idxlast = -1;

// Declaração global
Double_t phi_last_quark, eta_last_quark, pt_last_quark;

Int_t findLastQuark(TClonesArray *particulas, Int_t index = -1);

void juntandoBBCC()
{
    std::ofstream file("pesadosBar_Tagging.csv");
    file << "PtJet,nConstituent,Angularidade,SigmaKT,MaxRho,N_Vertices_S,FracLeadPt,Label" << std::endl;


    std::vector<std::string> arquivos = {"BTagging.root", "CTagging.root"};

    TH1F *histPtJato = new TH1F("histPtJato","Pt jato bottom; p_{T} [GeV]; Eventos", 100, 0, 80);
    TH1F *histPtJatoNao = new TH1F("histPtJatoNao","Pt jato nao-bottom; p_{T} [GeV]; Eventos", 100, 0, 80);

    TH1F *histNConstitBottom = new TH1F("histNConstitBottom","N Constituintes Bottom; N_{const}; Eventos", 40, 0, 20);
    TH1F *histNConstitNonBottom = new TH1F("histNConstitNonBottom","N Constituintes Non-Bottom; N_{const}; Eventos", 40, 0, 20);

    TH1F *histAngBottom = new TH1F("histAngBottom","Ângulo Bottom; Ângulo [rad]; Eventos", 100, 0, 4);
    TH1F *histAngNonBottom = new TH1F("histAngNonBottom","Ângulo Non-Bottom; Ângulo [rad]; Eventos", 100, 0, 4);

    TH1F *histSigmaKTBottom = new TH1F("histSigmaKTBottom","#sigma_{k_{T}} Bottom; #sigma_{k_{T}} [GeV]; Eventos", 100, 0, 40);
    TH1F *histSigmaKTNonBottom = new TH1F("histSigmaKTNonBottom","#sigma_{k_{T}} Non-Bottom; #sigma_{k_{T}} [GeV]; Eventos", 100, 0, 40);

    TH1F *histMaxRhoBottom = new TH1F("histMaxRhoBottom","Rho Max Bottom; #rho_{max}; Eventos", 100, 0, 2.1);
    TH1F *histMaxRhoNonBottom = new TH1F("histMaxRhoNonBottom","Rho Max Non-Bottom; #rho_{max}; Eventos", 100, 0, 2.1);

    TH1F *histNSecVBottom = new TH1F("histNSecVBottom","N Secundarias Verticais Bottom; N_{secV}; Eventos", 15, 0, 15);
    TH1F *histNSecVNonBottom = new TH1F("histNSecVNonBottom","N Secundarias Verticais Non-Bottom; N_{secV}; Eventos", 15, 0, 15);

    TH1F *histFracLeadPtBottom = new TH1F("histFracLeadPtBottom", "Frac Lead Pt Bottom; p_{T}^{lead}/p_{T}^{jet}; Eventos", 100, 1, 1);
    TH1F *histFracLeadPtNonBottom = new TH1F("histFracLeadPtNonBottom","Frac Lead Pt Non-Bottom; p_{T}^{lead}/p_{T}^{jet}; Eventos", 100, 1, 1);



    for (auto &nomeArquivo : arquivos) {
    TFile *rootfile = TFile::Open(nomeArquivo.c_str(), "READ");

    if (!rootfile || !rootfile->IsOpen()) {
        std::cerr << "Erro: Não consegui abrir " << nomeArquivo << std::endl;
        continue;
    }

    TTree *ttree = dynamic_cast<TTree *>(rootfile->Get("tree"));
    if (!ttree) {
        std::cerr << "Erro: Não achei TTree em " << nomeArquivo << std::endl;
        rootfile->Close();
        continue;
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

        // Lista de PDG codes para os hádrons de c e b que serão aceitos
        const int nBottomHadrons = 62;
    int BottomHadrons[nBottomHadrons] = {511, 521, 513, 523, 531, 533, 541, 543, 5122, 5112, 5212, 5222, 5132, 5232, 5332, 5142, 5242, 5342, 5442, 5512, 5522, 5532, 5542, 5554, 5152, 5252, 5352, 5162, 5262, 5362, 5454,
                                            411, 421, 413, 423, 415, 425, 431, 433, 435, 4112, 4212, 4222, 4122, 4132, 4232, 4332, 4312, 4322, 4114, 4214, 4224, 4134, 4234, 4334, 4412, 4422, 4432, 4414, 4424, 4434, 4444};


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

                           
                            for (int i = 0; i < nBottomHadrons; i++)
                            {

                                if (absMotherPdg == BottomHadrons[i])
                                {
                                    MyPart *mp = static_cast<MyPart *>(particulas->At(ip));

                                    // pega o pdg diretamente de mp
                                    int pdg = mp->fPdg;

                                    // Int_t idxQuarkAntesHadron = -1;
                                
                                    Int_t idxLastCharmQuark = -1;
                                    if (traceMotherForBottomHadron(particulas, motherIndex, idxQuarkAntesHadron))
                                    {

                                        // Para utilizar na construcao do match
                                        idxlast = idxQuarkAntesHadron;
                                        jet_has_bottom = true;
                                        break; 
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
                            // Double_t PDG = quarkBottom->fPdg;

                            // std::cout << "pdg: " << PDG << std::endl;

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
                                
                                histPtJato->Fill(jetPt);
                                histNConstitBottom->Fill(nConstituents);
                                histAngBottom->Fill(Angularidade);
                                histSigmaKTBottom->Fill(DesvioSigmaKt);
                                histMaxRhoBottom->Fill(MaxRho);
                                histNSecVBottom->Fill(NSecV);
                                histFracLeadPtBottom->Fill(fracLeadPt);
                            }
                        }

                            
                    }
                    else
                    {
                        // Caso nao seja um jato de charm
                        Int_t Label = 0;
                        file << jetPt << "," << nConstituents << "," << Angularidade << "," << DesvioSigmaKt << "," << MaxRho << "," << NSecV << "," << fracLeadPt << "," << Label << std::endl;

                        histPtJatoNao->Fill(jetPt);
                        histNConstitNonBottom->Fill(nConstituents);
                        histAngNonBottom->Fill(Angularidade);
                        histSigmaKTNonBottom->Fill(DesvioSigmaKt);
                        histNSecVNonBottom->Fill(NSecV);
                        histMaxRhoNonBottom->Fill(MaxRho); 
                        histFracLeadPtNonBottom->Fill(fracLeadPt); 
                    }
                }

            } // Fim do loop dos jatos

        } // Fim do loop dos eventos

        

    rootfile->Close();
    }

    file.close();

    std::vector<std::pair<TH1F*, TH1F*>> hists = {
    {histPtJato,        histPtJatoNao},
    {histNConstitBottom, histNConstitNonBottom},
    {histAngBottom,      histAngNonBottom},
    {histSigmaKTBottom,  histSigmaKTNonBottom},
    {histMaxRhoBottom,   histMaxRhoNonBottom},
    {histNSecVBottom,    histNSecVNonBottom},
    {histFracLeadPtBottom, histFracLeadPtNonBottom}
    };

    int canvasIndex = 0;
    for (auto &pair : hists) {
        TH1F *hBottom    = pair.first;
        TH1F *hNonBottom = pair.second;

      
        if (hBottom->Integral() != 0)
            hBottom->Scale(1.0 / hBottom->Integral());
        if (hNonBottom->Integral() != 0)
            hNonBottom->Scale(1.0 / hNonBottom->Integral());

       
        TString cname = Form("canvas_%d", canvasIndex++);
        TCanvas *c = new TCanvas(cname, hBottom->GetTitle(), 800, 600);

      
        hBottom->SetLineColor(kBlue);
        hBottom->SetLineWidth(2);

        hNonBottom->SetLineColor(kRed);
        hNonBottom->SetLineWidth(2);
        hNonBottom->SetLineStyle(2);

       
        hBottom->SetYTitle("Normalizado");

      
        hBottom->Draw("HIST");
        hNonBottom->Draw("HIST SAME");

        
        auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
        legend->AddEntry(hBottom,    "Jatos com bottom", "l");
        legend->AddEntry(hNonBottom, "Jatos sem bottom", "l");
        legend->Draw();
    }

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

    if ((idxLastCharmQuark == -1) && ((TMath::Abs(mom->fPdg) == 5) || (TMath::Abs(mom->fPdg) == 4)))
    {
        idxLastCharmQuark = momIdx;
    }

    // Parou num quark charm?
    if ((momIdx == 4 || momIdx == 5) && ((TMath::Abs(mom->fPdg) == 5) || (TMath::Abs(mom->fPdg) == 4)))
    {

        return true;
    }

    // Senão, sobe recursivamente
    return traceMotherForBottomHadron(particulas, momIdx, idxLastCharmQuark);
}
