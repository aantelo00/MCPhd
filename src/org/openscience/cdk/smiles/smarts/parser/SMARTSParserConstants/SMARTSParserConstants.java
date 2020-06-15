package org.openscience.cdk.smiles.smarts.parser.SMARTSParserConstants;

//FIXMED CLASE VACIA
public interface SMARTSParserConstants {
	 int EOF = 0;
	  int _WS = 1;
	  int L_AND = 2;
	  int H_AND = 3;
	  int OR = 4;
	  int NOT = 5;
	  int S_BOND = 6;
	  int UP_S_BOND = 7;
	  int DN_S_BOND = 8;
	  int UP_OR_UNSPECIFIED_S_BOND = 9;
	  int DN_OR_UNSPECIFIED_S_BOND = 10;
	  int D_BOND = 11;
	  int T_BOND = 12;
	  int AR_BOND = 13;
	  int ANY_BOND = 14;
	  int R_BOND = 15;
	  int c = 16;
	  int n = 17;
	  int o = 18;
	  int s = 19;
	  int p = 20;
	  int as = 21;
	  int se = 22;
	  int H = 23;
	  int HE = 24;
	  int LI = 25;
	  int BE = 26;
	  int B = 27;
	  int C = 28;
	  int N = 29;
	  int O = 30;
	  int F = 31;
	  int NE = 32;
	  int NA = 33;
	  int MG = 34;
	  int AL = 35;
	  int SI = 36;
	  int P = 37;
	  int S = 38;
	  int CL = 39;
	  int AR = 40;
	  int K = 41;
	  int CA = 42;
	  int SC = 43;
	  int TI = 44;
	  int V = 45;
	  int CR = 46;
	  int MN = 47;
	  int FE = 48;
	  int CO = 49;
	  int NI = 50;
	  int CU = 51;
	  int ZN = 52;
	  int GA = 53;
	  int GE = 54;
	  int AS = 55;
	  int SE = 56;
	  int BR = 57;
	  int KR = 58;
	  int RB = 59;
	  int SR = 60;
	  int Y = 61;
	  int ZR = 62;
	  int NB = 63;
	  int MO = 64;
	  int TC = 65;
	  int RU = 66;
	  int RH = 67;
	  int PD = 68;
	  int AG = 69;
	  int CD = 70;
	  int IN = 71;
	  int SN = 72;
	  int SB = 73;
	  int TE = 74;
	  int I = 75;
	  int XE = 76;
	  int CS = 77;
	  int BA = 78;
	  int LA = 79;
	  int HF = 80;
	  int TA = 81;
	  int W = 82;
	  int RE = 83;
	  int OS = 84;
	  int IR = 85;
	  int PT = 86;
	  int AU = 87;
	  int HG = 88;
	  int TL = 89;
	  int PB = 90;
	  int BI = 91;
	  int PO = 92;
	  int AT = 93;
	  int RN = 94;
	  int FR = 95;
	  int RA = 96;
	  int AC = 97;
	  int TH = 98;
	  int PA = 99;
	  int U = 100;
	  int PU = 101;
	  int AM = 102;
	  int CM = 103;
	  int BK = 104;
	  int CF = 105;
	  int ES = 106;
	  int FM = 107;
	  int MD = 108;
	  int NO = 109;
	  int LR = 110;
	  int NP = 111;
	  int CE = 112;
	  int ND = 113;
	  int PM = 114;
	  int SM = 115;
	  int EU = 116;
	  int GD = 117;
	  int TB = 118;
	  int DY = 119;
	  int HO = 120;
	  int ER = 121;
	  int TM = 122;
	  int YB = 123;
	  int LU = 124;
	  int PR = 125;
	  int WILDCARD = 126;
	  int h = 127;
	  int a = 128;
	  int A = 129;
	  int D = 130;
	  int R = 131;
	  int r = 132;
	  int v = 133;
	  int X = 134;
	  int x = 135;
	  int G = 136;
	  int HX = 137;
	  int CARET = 138;
	  int DOLLAR = 139;
	  int L_PAREN = 140;
	  int R_PAREN = 141;
	  int L_BRACKET = 142;
	  int R_BRACKET = 143;
	  int Q_MARK = 144;
	  int DIGIT = 145;

	  int DEFAULT = 0;

	  String[] tokenImage = {
	    "<EOF>",
	    "<_WS>",
	    "\";\"",
	    "\"&\"",
	    "\",\"",
	    "\"!\"",
	    "\"-\"",
	    "\"/\"",
	    "\"\\\\\"",
	    "\"/?\"",
	    "\"\\\\?\"",
	    "\"=\"",
	    "\"#\"",
	    "\":\"",
	    "\"~\"",
	    "\"@\"",
	    "\"c\"",
	    "\"n\"",
	    "\"o\"",
	    "\"s\"",
	    "\"p\"",
	    "\"as\"",
	    "\"se\"",
	    "\"H\"",
	    "\"He\"",
	    "\"Li\"",
	    "\"Be\"",
	    "\"B\"",
	    "\"C\"",
	    "\"N\"",
	    "\"O\"",
	    "\"F\"",
	    "\"Ne\"",
	    "\"Na\"",
	    "\"Mg\"",
	    "\"Al\"",
	    "\"Si\"",
	    "\"P\"",
	    "\"S\"",
	    "\"Cl\"",
	    "\"Ar\"",
	    "\"K\"",
	    "\"Ca\"",
	    "\"Sc\"",
	    "\"Ti\"",
	    "\"V\"",
	    "\"Cr\"",
	    "\"Mn\"",
	    "\"Fe\"",
	    "\"Co\"",
	    "\"Ni\"",
	    "\"Cu\"",
	    "\"Zn\"",
	    "\"Ga\"",
	    "\"Ge\"",
	    "\"As\"",
	    "\"Se\"",
	    "\"Br\"",
	    "\"Kr\"",
	    "\"Rb\"",
	    "\"Sr\"",
	    "\"Y\"",
	    "\"Zr\"",
	    "\"Nb\"",
	    "\"Mo\"",
	    "\"Tc\"",
	    "\"Ru\"",
	    "\"Rh\"",
	    "\"Pd\"",
	    "\"Ag\"",
	    "\"Cd\"",
	    "\"In\"",
	    "\"Sn\"",
	    "\"Sb\"",
	    "\"Te\"",
	    "\"I\"",
	    "\"Xe\"",
	    "\"Cs\"",
	    "\"Ba\"",
	    "\"La\"",
	    "\"Hf\"",
	    "\"Ta\"",
	    "\"W\"",
	    "\"Re\"",
	    "\"Os\"",
	    "\"Ir\"",
	    "\"Pt\"",
	    "\"Au\"",
	    "\"Hg\"",
	    "\"Tl\"",
	    "\"Pb\"",
	    "\"Bi\"",
	    "\"Po\"",
	    "\"At\"",
	    "\"Rn\"",
	    "\"Fr\"",
	    "\"Ra\"",
	    "\"Ac\"",
	    "\"Th\"",
	    "\"Pa\"",
	    "\"U\"",
	    "\"Pu\"",
	    "\"Am\"",
	    "\"Cm\"",
	    "\"Bk\"",
	    "\"Cf\"",
	    "\"Es\"",
	    "\"Fm\"",
	    "\"Md\"",
	    "\"No\"",
	    "\"Lr\"",
	    "\"Np\"",
	    "\"Ce\"",
	    "\"Nd\"",
	    "\"Pm\"",
	    "\"Sm\"",
	    "\"Eu\"",
	    "\"Gd\"",
	    "\"Tb\"",
	    "\"Dy\"",
	    "\"Ho\"",
	    "\"Er\"",
	    "\"Tm\"",
	    "\"Yb\"",
	    "\"Lu\"",
	    "\"Pr\"",
	    "\"*\"",
	    "\"h\"",
	    "\"a\"",
	    "\"A\"",
	    "\"D\"",
	    "\"R\"",
	    "\"r\"",
	    "\"v\"",
	    "\"X\"",
	    "\"x\"",
	    "\"G\"",
	    "\"#X\"",
	    "\"^\"",
	    "\"$\"",
	    "\"(\"",
	    "\")\"",
	    "\"[\"",
	    "\"]\"",
	    "\"?\"",
	    "<DIGIT>",
	    "\">>\"",
	    "\">\"",
	    "\".\"",
	    "\"+\"",
	    "\"--\"",
	    "\"---\"",
	    "\"----\"",
	    "\"-----\"",
	    "\"------\"",
	    "\"-------\"",
	    "\"--------\"",
	    "\"++\"",
	    "\"+++\"",
	    "\"++++\"",
	    "\"+++++\"",
	    "\"++++++\"",
	    "\"+++++++\"",
	    "\"++++++++\"",
	    "\"@@\"",
	  };

	}
