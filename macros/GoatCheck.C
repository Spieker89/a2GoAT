

void GoatCheck(const char* dirname, const char* prefix)
{
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) 
   {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      TH2D* ChiSqDist = 0;
      TH2D* TOF = 0;
      
      printf("rm ");
      
      while ((file=(TSystemFile*)next())) 
      {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.BeginsWith(prefix) && fname.EndsWith(".root")) 
         {
			char name[1024];
			sprintf(name, "%s/%s", dirname, file->GetName());
			//cout << name << endl;
			TFile *f1 = TFile::Open(name);
			if(!f1) continue;
			if(!f1->IsOpen()) continue;
			if(!f1->Get("protons"))
			{
				printf("%s ", name);
				continue;
			}
			Int_t e = ((TTree*)f1->Get("protons"))->GetEntries();
			if(e != ((TTree*)f1->Get("etaPrimes"))->GetEntries())
				printf("********Error etap\n");
			if(e != ((TTree*)f1->Get("eventParameters"))->GetEntries())
				printf("********Error eventParameters\n");
			if(e != ((TTree*)f1->Get("photons"))->GetEntries())
				printf("********Error photons\n");
			if(e != ((TTree*)f1->Get("tracks"))->GetEntries())
				printf("********Error tracks\n");
			if(e != ((TTree*)f1->Get("tagger"))->GetEntries())
				printf("********Error tagger\n");
		}
      }
   }
}
