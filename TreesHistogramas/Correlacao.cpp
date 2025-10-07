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

Int_t idxlast = -1;

// Declaração global
Double_t phi_last_quark, eta_last_quark, pt_last_quark;

Int_t findLastQuark(TClonesArray *particulas, Int_t index = -1);

void Correlacao()
{
    TFile *rootfile = TFile::Open("HardQCDO.root", "READ");

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

    TH2F *histPtNConstit       = new TH2F("histPtNConstit", "p_{T} do jato vs N^{o} Constituintes", 100, 0, 40, 40, 0, 20);
    TH2F *histPtSigmaKT        = new TH2F("histPtSigmaKT", "p_{T} do jato vs #sigma_{k_{T}}", 100, 0, 40, 100, 0, 10);
    TH2F *histPtFracLeadPt     = new TH2F("histPtFracLeadPt", "p_{T} do jato vs Fracao Lead", 100, 0, 40, 100, 0, 1);
    TH2F *histPtNSecV          = new TH2F("histPtNSecV", "p_{T} do jato vs N^{o} Vertices Secundarios", 100, 0, 40, 15, 0, 15);

    TH2F *histNConstitSigmaKT  = new TH2F("histNConstitSigmaKT", "N^{o} Constituintes vs #sigma_{k_{T}}", 40, 0, 20, 100, 0, 10);
    TH2F *histNConstitFracLead = new TH2F("histNConstitFracLead", "N^{o} Constituintes vs Fracao Lead", 40, 0, 20, 100, 0, 1);
    TH2F *histNConstitNSecV    = new TH2F("histNConstitNSecV", "N^{o} Constituintes vs N^{o} Vertices Secundarios", 40, 0, 20, 15, 0, 15);

    TH2F *histSigmaKTFracLead  = new TH2F("histSigmaKTFracLead", "#sigma_{k_{T}} vs Fracao Lead", 100, 0, 10, 100, 0, 1);
    TH2F *histSigmaKTNSecV     = new TH2F("histSigmaKTNSecV", "#sigma_{k_{T}} vs N^{o} Vertices Secundarios", 100, 0, 10, 15, 0, 15);

    TH2F *histFracLeadNSecV    = new TH2F("histFracLeadNSecV", "Fracao Lead vs N^{o} Vertices Secundarios", 100, 0, 1, 15, 0, 15);


    
    TClonesArray *particulas = new TClonesArray("MyPart");
    TClonesArray *jatos = new TClonesArray("MyJet");
    TClonesArray *constituintes = new TClonesArray("MyConstituent");

    
    ttree->SetBranchAddress("particulas", &particulas);
    ttree->SetBranchAddress("jatos", &jatos);
    ttree->SetBranchAddress("constituintes", &constituintes);

    Long64_t nev = ttree->GetEntries();


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


                //---------------------------------------------------------------------------------------------------------
                // Variaveis derivadas utilizadas no aprendizado de maquinas
                //---------------------------------------------------------------------------------------------------------

                //---------------------------------------------------------------------------------------------------------
                // Cálculo de Rho (distância radial no plano XY)
                Double_t Rho = TMath::Sqrt(pow(vx, 2) + pow(vy, 2));  

                if ((Rho > 0.01 ) && (Rho <= 1))
                {
                    NSecV++;
                    if (Rho > MaxRho){
                        
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

                
                Double_t pt = jetPt;
                Double_t nConstit = nConstituents;
                Double_t sigma = DesvioSigmaKt;
                Double_t fracLead = fracLeadPt;
                Double_t nSecV = NSecV;

               
                histPtNConstit->Fill(jetPt, nConstituents);
                histPtSigmaKT->Fill(jetPt, DesvioSigmaKt);
                histPtFracLeadPt->Fill(jetPt, fracLeadPt);
                histPtNSecV->Fill(jetPt, NSecV);

                histNConstitSigmaKT->Fill(nConstituents, DesvioSigmaKt);
                histNConstitFracLead->Fill(nConstituents, fracLeadPt);
                histNConstitNSecV->Fill(nConstituents, NSecV);

                histSigmaKTFracLead->Fill(DesvioSigmaKt, fracLeadPt);
                histSigmaKTNSecV->Fill(DesvioSigmaKt, NSecV);

                histFracLeadNSecV->Fill(fracLeadPt, NSecV);
                
            }

        } // Fim do loop dos jatos

    } // Fim do loop dos eventos


gStyle->SetOptStat(1110);         
gStyle->SetStatFontSize(0.06);   
gStyle->SetStatX(0.88);           
gStyle->SetStatY(0.92);          
gStyle->SetStatW(0.28);           
gStyle->SetStatH(0.22);           




TCanvas *cAll = new TCanvas("cAll", "Todas as correlacoes", 2400, 1600);
cAll->Divide(4,3); 


auto configHist = [](TH2F* h, const char* xtitle, const char* ytitle){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

  
    h->GetXaxis()->SetTitleSize(0.07);
    h->GetYaxis()->SetTitleSize(0.07);

    
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);

  
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetZaxis()->SetLabelSize(0.055); 

    
    h->SetStats(1);
};


auto configPad = [](){
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.08);
};


cAll->cd(1);  configPad();
configHist(histPtNConstit, "p_{T} [GeV]", "N_{const}");
histPtNConstit->Draw("COLZ");

cAll->cd(2);  configPad();
configHist(histPtSigmaKT, "p_{T} [GeV]", "#sigma_{kT} [GeV]");
histPtSigmaKT->Draw("COLZ");

cAll->cd(3);  configPad();
configHist(histPtFracLeadPt, "p_{T} [GeV]", "p_{T}^{lead}/p_{T}^{jet}");
histPtFracLeadPt->Draw("COLZ");

cAll->cd(4);  configPad();
configHist(histPtNSecV, "p_{T} [GeV]", "N_{secV}");
histPtNSecV->Draw("COLZ");


cAll->cd(5);  configPad();
configHist(histNConstitSigmaKT, "N_{const}", "#sigma_{kT} [GeV]");
histNConstitSigmaKT->Draw("COLZ");

cAll->cd(6);  configPad();
configHist(histNConstitFracLead, "N_{const}", "p_{T}^{lead}/p_{T}^{jet}");
histNConstitFracLead->Draw("COLZ");

cAll->cd(7);  configPad();
configHist(histNConstitNSecV, "N_{const}", "N_{secV}");
histNConstitNSecV->Draw("COLZ");

cAll->cd(8);  configPad();
configHist(histSigmaKTFracLead, "#sigma_{kT} [GeV]", "p_{T}^{lead}/p_{T}^{jet}");
histSigmaKTFracLead->Draw("COLZ");


cAll->cd(9);  configPad();
configHist(histSigmaKTNSecV, "#sigma_{kT} [GeV]", "N_{secV}");
histSigmaKTNSecV->Draw("COLZ");

cAll->cd(10); configPad();
configHist(histFracLeadNSecV, "p_{T}^{lead}/p_{T}^{jet}", "N_{secV}");
histFracLeadNSecV->Draw("COLZ");


   
}

