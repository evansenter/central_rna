/*

  prototypes for energy_par.c

*/
extern double triplet37[256][16];
extern double stack37[16][16]; // account for AU,AC,AG,GC,GU,GA,CG,CU,CA,UA,UG,UC
extern double hairpin37[31];
extern double hairpin3_37[16][64];
extern double hairpin4_37[16][256];
extern double hairpin6_37[256][256];
extern double bulge37[31];
extern double internal37[31];
extern double internal11_37[256][16];
extern double internal12_37[256][64];
extern double internal22_37[256][256];
extern double dangle37[16][8]; // i - AU,AC,AG,GC,GU,GA,CG,CU,CA,UA,UG,UC 
                               // j - 0A,A0,0U,U0,0C,C0,0G,G0
extern double TMM37[16][16]; // have to include GG, AA, CC, UU as well

extern double AUGUinternal37;
extern double intern_init37; // intermolecular initiation  
extern double C_bulge37; // special C bulge penalty
extern double GUAU_penalty37; // penalty for ending helix with AU or GU
extern double GGmismatch37;
extern double GGGhairpin37;
extern double UUGAmismatch37;
extern double GUclosure37;
extern double C3loop37;
extern double CloopA37;
extern double CloopB37;
extern double Initiation37;
extern double assymetry37;
extern double MultiloopStrain37;
extern double MultiloopA37;
extern double MultiloopB37;
extern double MultiloopC37;
extern double MultiloopAssym37;
extern double DisconPenalty37;

extern double triplet[256][16];
extern double stack[16][16]; // account for AU,AC,AG,GC,GU,GA,CG,CU,CA,UA,UG,UC
extern double hairpin[31];
extern double hairpin3[16][64];
extern double hairpin4[16][256];
extern double hairpin6[256][256];
extern double bulge[31];
extern double internal[31];
extern double internal11[256][16];
extern double internal12[256][64];
extern double internal22[256][256];
extern double dangle[16][8]; // i - AU,AC,AG,GC,GU,GA,CG,CU,CA,UA,UG,UC 
                               // j - 0A,A0,0U,U0,0C,C0,0G,G0
extern double TMM[16][16]; // have to include GG, AA, CC, UU as well

extern double intern_init; // intermolecular initiation  
extern double C_bulge; // special C bulge penalty
extern double GUAU_penalty; // penalty for ending helix with AU or GU
extern double GGmismatch;
extern double GGGhairpin;
extern double assymetry;
extern double UUGAmismatch;
extern double GUclosure;
extern double AUGUinternal;
extern double C3loop;
extern double CloopA;
extern double CloopB;
extern double Initiation;
extern double MultiloopStrain;
extern double MultiloopA;
extern double MultiloopB;
extern double MultiloopC;
extern double MultiloopAssym;
extern double DisconPenalty;
