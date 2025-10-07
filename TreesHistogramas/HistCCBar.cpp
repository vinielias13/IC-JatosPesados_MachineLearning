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

bool traceMothersForCharmHadron(TClonesArray *particles, Int_t index, Int_t &idxLastCharmQuark);

Int_t idxlast = -1;

// Declaração global
Double_t phi_last_quark, eta_last_quark, pt_last_quark;

Int_t findLastQuark(TClonesArray *particulas, Int_t index = -1);

void HistCCBar()
{

    TFile *rootfile = TFile::Open("CTagging.root", "READ");

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

    TH1F *histPtJato = new TH1F("histPtJato","Distribuicao do p_{T} do jato charm; p_{T} [GeV]; Eventos", 100, 0, 80);
    TH1F *histPtJatoNao = new TH1F("histPtJatoNao","Distribuicao do p_{T} do jato Nao-Charm; p_{T} [GeV]; Eventos", 100, 0, 80);

    TH1F *histNConstitCharm = new TH1F("histNConstitCharm", "N^{o}  de Constituintes Charm; N_{const}; Eventos", 40, 0, 20);
    TH1F *histNConstitNonCharm = new TH1F("histNConstitNonCharm","N^{o}  de Constituintes Nao-Charm; N_{const}; Eventos", 40, 0, 20);

    TH1F *histAngCharm = new TH1F("histAngCharm", "Angulo Charm; Angulo (rad); Eventos", 100, 0, 4);
    TH1F *histAngNonCharm = new TH1F("histAngNonCharm", "Angulo Non-Charm; Angulo (rad); Eventos", 100, 0, 4);

    TH1F *histSigmaKTCharm = new TH1F("histSigmaKTCharm",    "#sigma_{k_{T}} Charm; #sigma_{k_{T}} (GeV); Eventos", 100, 0, 10);
    TH1F *histSigmaKTNonCharm = new TH1F("histSigmaKTNonCharm", "#sigma_{k_{T}} Nao-Charm; #sigma_{k_{T}} (GeV); Eventos", 100, 0, 10);

    TH1F *histMaxRhoCharm = new TH1F("histMaxRhoCharm", "Rho Max dos Jatos Charm e Nao-Charm; mm", 100, 1, 1);
    TH1F *histMaxRhoNonCharm = new TH1F("histMaxRho", "Rho Max  dos Jatos Charm e Nao-Charm; mm", 100, 1, 1);

    TH1F *histNSecVCharm = new TH1F("histNSecVCharm", "N^{o}  Vertices Secundarios Charm; N_{secV}; Eventos", 15, 0, 15);
    TH1F *histNSecVNonCharm = new TH1F("histNSecVNonCharm", "N^{o}  Vertices Secundarios Nao-Charm; N_{secV}; Eventos", 15, 0, 15);

    TH1F *histFracLeadPtCharm    = new TH1F("histFracLeadPtCharm","Fracao Lead Pt Charm; p_{T}^{lead}/p_{T}^{jet}; Eventos", 100, 1, 1);
    TH1F *histFracLeadPtNonCharm = new TH1F("histFracLeadPtNonCharm","Fracao Lead Pt Non-Charm; p_{T}^{lead}/p_{T}^{jet}; Eventos", 100, 1, 1);

    // TH2F *histPtNConstit       = new TH2F("histPtNConstit", "p_{T} do jato vs N^{o} Constituintes (Charm)", 100, 0, 40, 40, 0, 20);
    // TH2F *histPtSigmaKT        = new TH2F("histPtSigmaKT", "p_{T} do jato vs #sigma_{k_{T}} (Charm)", 100, 0, 40, 100, 0, 10);
    // TH2F *histPtFracLeadPt     = new TH2F("histPtFracLeadPt", "p_{T} do jato vs Fracao Lead (Charm)", 100, 0, 40, 100, 0, 1);
    // TH2F *histPtNSecV          = new TH2F("histPtNSecV", "p_{T} do jato vs N^{o} Vertices Secundarios (Charm)", 100, 0, 40, 15, 0, 15);

    // TH2F *histNConstitSigmaKT  = new TH2F("histNConstitSigmaKT", "N^{o} Constituintes vs #sigma_{k_{T}} (Charm)", 40, 0, 20, 100, 0, 10);
    // TH2F *histNConstitFracLead = new TH2F("histNConstitFracLead", "N^{o} Constituintes vs Fracao Lead (Charm)", 40, 0, 20, 100, 0, 1);
    // TH2F *histNConstitNSecV    = new TH2F("histNConstitNSecV", "N^{o} Constituintes vs N^{o} Vertices Secundarios (Charm)", 40, 0, 20, 15, 0, 15);

    // TH2F *histSigmaKTFracLead  = new TH2F("histSigmaKTFracLead", "#sigma_{k_{T}} vs Fracao Lead (Charm)", 100, 0, 10, 100, 0, 1);
    // TH2F *histSigmaKTNSecV     = new TH2F("histSigmaKTNSecV", "#sigma_{k_{T}} vs N^{o} Vertices Secundarios (Charm)", 100, 0, 10, 15, 0, 15);

    // TH2F *histFracLeadNSecV    = new TH2F("histFracLeadNSecV", "Fracao Lead vs N^{o} Vertices Secundarios (Charm)", 100, 0, 1, 15, 0, 15);

    // O Array com as variaveis de interesse
    TClonesArray *particulas = new TClonesArray("MyPart");
    TClonesArray *jatos = new TClonesArray("MyJet");
    TClonesArray *constituintes = new TClonesArray("MyConstituent");

    // Coluna com as variaveis
    ttree->SetBranchAddress("particulas", &particulas);
    ttree->SetBranchAddress("jatos", &jatos);
    ttree->SetBranchAddress("constituintes", &constituintes);

    Long64_t nev = ttree->GetEntries();

    // Lista de PDG codes para os hádrons de c que serão aceitos
    const int nCharmHadrons = 9;
    int charmHadrons[nCharmHadrons] = {411, 421, 413, 423, 415, 425, 431, 433, 435};

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


            // Boleano que se torna verdadeiro se o jato for proveninte de quark charm
            bool jet_has_charm = false;

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
                // Int_t ist = constituent->status;

                if (!jet_has_charm)
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

                        // Se achamos um hádron de charm na árvore genealógica, seguimos até o quark charm original
                        for (int i = 0; i < nCharmHadrons; i++)
                        {

                            // Se alguma mae estiver na lista dos hadrons de charm entra no If
                            if (absMotherPdg == charmHadrons[i])
                            {
                                MyPart *mp = static_cast<MyPart *>(particulas->At(ip));

                                // pega o pdg diretamente de mp
                                int pdg = mp->fPdg;

                                // Int_t idxQuarkAntesHadron = -1;
                                // Se a funcao retornar True sabeberemos que o jato e de charm
                                Int_t idxLastCharmQuark = -1;
                                if (traceMothersForCharmHadron(particulas, motherIndex, idxQuarkAntesHadron))
                                {

                                    // Para utilizar na construcao do match
                                    idxlast = idxQuarkAntesHadron;
                                    jet_has_charm = true;
                                    break; // Se econtrarmos pelo menos um hadron de charm saimos do loop
                                }

                            }

                        }

                        if (jet_has_charm)
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

                    if ((Rho <= 1) && (Rho > MaxRho))
                    {
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

            // Double_t DesvioSigmaKt = sigmaKT / nConstituents;

            //---------------------------------------------------------------------------------------------------------

            if ((jetPt > 5))
            {

                if (jet_has_charm)
                {   
                    if (idxlast >= 0){

                        MyPart *quarkCharm = static_cast<MyPart *>(particulas->At(idxlast));

                        double_t ptQuark = quarkCharm->fPt;
                        Double_t PhiQuark = quarkCharm->fPhi;
                        Double_t EtaQuark = quarkCharm->fEta;

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
                            
                            histPtJato->Fill(jetPt);            
                            histNConstitCharm->Fill(nConstituents);
                            histAngCharm->Fill(Angularidade);            
                            histSigmaKTCharm->Fill(DesvioSigmaKt); 
                            histMaxRhoCharm->Fill(MaxRho); 
                            histNSecVCharm->Fill(NSecV);         
                            histFracLeadPtCharm->Fill(fracLeadPt);

                            
                            // Double_t pt = jetPt;
                            // Double_t nConstit = nConstituents;
                            // Double_t sigma = DesvioSigmaKt;
                            // Double_t fracLead = fracLeadPt;
                            // Double_t nSecV = NSecV;

                            // // Preenchendo os histogramas 2D de correlação
                            // histPtNConstit->Fill(jetPt, nConstituents);
                            // histPtSigmaKT->Fill(jetPt, DesvioSigmaKt);
                            // histPtFracLeadPt->Fill(jetPt, fracLeadPt);
                            // histPtNSecV->Fill(jetPt, NSecV);

                            // histNConstitSigmaKT->Fill(nConstituents, DesvioSigmaKt);
                            // histNConstitFracLead->Fill(nConstituents, fracLeadPt);
                            // histNConstitNSecV->Fill(nConstituents, NSecV);

                            // histSigmaKTFracLead->Fill(DesvioSigmaKt, fracLeadPt);
                            // histSigmaKTNSecV->Fill(DesvioSigmaKt, NSecV);

                            // histFracLeadNSecV->Fill(fracLeadPt, NSecV);
                            
                        }
                    }
                        
                }
                else
                {
                    // Caso nao seja um jato de charm
                    histPtJatoNao->Fill(jetPt);            
                    histNConstitNonCharm->Fill(nConstituents);
                    histAngNonCharm->Fill(Angularidade);            
                    histSigmaKTNonCharm->Fill(DesvioSigmaKt);
                    histMaxRhoNonCharm->Fill(MaxRho);       
                    histNSecVNonCharm->Fill(NSecV);   
                    histFracLeadPtNonCharm->Fill(fracLeadPt);   
           
                
                }
            }

        } // Fim do loop dos jatos

    } // Fim do loop dos eventos

    TCanvas *cAll = new TCanvas("cAll", "Distribuicoes - Charm vs Nao-Charm", 1600, 1000);
cAll->Divide(3, 2); 


std::vector<std::pair<TH1F*, TH1F*>> selectedHists = {
    {histPtJato,        histPtJatoNao},
    {histNConstitCharm, histNConstitNonCharm},
    {histSigmaKTCharm,  histSigmaKTNonCharm},
    {histFracLeadPtCharm, histFracLeadPtNonCharm},
    {histNSecVCharm,    histNSecVNonCharm}
    
};


gStyle->SetOptStat(1110); 
gStyle->SetStatColor(0);  

int padIndex = 1;
for (auto &pair : selectedHists) {
    TH1F *hCharm    = pair.first;
    TH1F *hNonCharm = pair.second;

    
    if (hCharm->Integral() != 0)    hCharm->Scale(1.0 / hCharm->Integral());
    if (hNonCharm->Integral() != 0) hNonCharm->Scale(1.0 / hNonCharm->Integral());

    
    cAll->cd(padIndex);

    hCharm->SetLineColor(kBlue);
    hCharm->SetLineWidth(2);

    hNonCharm->SetLineColor(kRed);
    hNonCharm->SetLineWidth(2);
    hNonCharm->SetLineStyle(2);

    hCharm->SetYTitle("Normalizado");

    
    double maxY = std::max(hCharm->GetMaximum(), hNonCharm->GetMaximum());
    hCharm->SetMaximum(1.2 * maxY); 

    // Desenha Charm e Non-Charm
    hCharm->Draw("HIST");
    hNonCharm->Draw("HIST SAMES");

    gPad->Update(); 

    // Caixa Charm
    TPaveStats *stCharm = (TPaveStats*) hCharm->GetListOfFunctions()->FindObject("stats");
    if (stCharm) {
        stCharm->SetTextColor(kBlue);
        stCharm->SetLineColor(kBlue);
        stCharm->SetX1NDC(0.15);
        stCharm->SetX2NDC(0.45);
        stCharm->SetY1NDC(0.75);
        stCharm->SetY2NDC(0.9);
    }

    // Caixa Não-Charm
    TPaveStats *stNonCharm = (TPaveStats*) hNonCharm->GetListOfFunctions()->FindObject("stats");
    if (stNonCharm) {
        stNonCharm->SetTextColor(kRed);
        stNonCharm->SetLineColor(kRed);
        stNonCharm->SetX1NDC(0.55);
        stNonCharm->SetX2NDC(0.85);
        stNonCharm->SetY1NDC(0.75);
        stNonCharm->SetY2NDC(0.9);
    }

    // Legenda
    auto legend = new TLegend(0.55, 0.6, 0.88, 0.73);
    legend->AddEntry(hCharm,    "Charm", "l");
    legend->AddEntry(hNonCharm, "Nao-Charm", "l");
    legend->Draw();

    gPad->Modified();
    gPad->Update();

    padIndex++;
}

    // //------------------------------------------------------- Para Correlacao ------------------------------------------------------------
    // gStyle->SetOptStat(1110);         
    // gStyle->SetStatFontSize(0.06);   
    // gStyle->SetStatX(0.88);           
    // gStyle->SetStatY(0.92);           
    // gStyle->SetStatW(0.28);           
    // gStyle->SetStatH(0.22);           
    // //gStyle->SetStatStyle(0);          


    // TCanvas *cAll = new TCanvas("cAll", "Todas as correlacoes", 2400, 1600);
    // cAll->Divide(4,3); 

    
    // auto configHist = [](TH2F* h, const char* xtitle, const char* ytitle){
    //     h->GetXaxis()->SetTitle(xtitle);
    //     h->GetYaxis()->SetTitle(ytitle);

    
    //     h->GetXaxis()->SetTitleSize(0.07);
    //     h->GetYaxis()->SetTitleSize(0.07);

    
    //     h->GetXaxis()->SetLabelSize(0.05);
    //     h->GetYaxis()->SetLabelSize(0.05);

   
    //     h->GetXaxis()->SetTitleOffset(1.2);
    //     h->GetYaxis()->SetTitleOffset(1.2);

    //     h->GetZaxis()->SetLabelSize(0.04); 

    //     // Caixa de estatística ligada
    //     h->SetStats(1);
    // };

    // auto configPad = [](){
    //     gPad->SetLeftMargin(0.18);
    //     gPad->SetBottomMargin(0.18);
    //     gPad->SetRightMargin(0.12);
    //     gPad->SetTopMargin(0.08);
    // };

    // cAll->cd(1);  configPad();
    // configHist(histPtNConstit, "p_{T} [GeV]", "N_{const}");
    // histPtNConstit->Draw("COLZ");

    // cAll->cd(2);  configPad();
    // configHist(histPtSigmaKT, "p_{T} [GeV]", "#sigma_{kT} [GeV]");
    // histPtSigmaKT->Draw("COLZ");

    // cAll->cd(3);  configPad();
    // configHist(histPtFracLeadPt, "p_{T} [GeV]", "p_{T}^{lead}/p_{T}^{jet}");
    // histPtFracLeadPt->Draw("COLZ");

    // cAll->cd(4);  configPad();
    // configHist(histPtNSecV, "p_{T} [GeV]", "N_{secV}");
    // histPtNSecV->Draw("COLZ");


    // cAll->cd(5);  configPad();
    // configHist(histNConstitSigmaKT, "N_{const}", "#sigma_{kT} [GeV]");
    // histNConstitSigmaKT->Draw("COLZ");

    // cAll->cd(6);  configPad();
    // configHist(histNConstitFracLead, "N_{const}", "p_{T}^{lead}/p_{T}^{jet}");
    // histNConstitFracLead->Draw("COLZ");

    // cAll->cd(7);  configPad();
    // configHist(histNConstitNSecV, "N_{const}", "N_{secV}");
    // histNConstitNSecV->Draw("COLZ");

    // cAll->cd(8);  configPad();
    // configHist(histSigmaKTFracLead, "#sigma_{kT} [GeV]", "p_{T}^{lead}/p_{T}^{jet}");
    // histSigmaKTFracLead->Draw("COLZ");


    // cAll->cd(9);  configPad();
    // configHist(histSigmaKTNSecV, "#sigma_{kT} [GeV]", "N_{secV}");
    // histSigmaKTNSecV->Draw("COLZ");

    // cAll->cd(10); configPad();
    // configHist(histFracLeadNSecV, "p_{T}^{lead}/p_{T}^{jet}", "N_{secV}");
    // histFracLeadNSecV->Draw("COLZ");


   
}

bool traceMothersForCharmHadron(TClonesArray *particulas, Int_t idx, Int_t &idxLastCharmQuark)
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

    if ((idxLastCharmQuark == -1) && (TMath::Abs(mom->fPdg) == 4))
    {
        idxLastCharmQuark = momIdx;
    }

    // Parou num quark charm?
    if ((momIdx == 4 || momIdx == 5) && (TMath::Abs(mom->fPdg) == 4))
    {

        return true;
    }

    // Senão, sobe recursivamente
    return traceMothersForCharmHadron(particulas, momIdx, idxLastCharmQuark);
}

