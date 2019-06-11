/* ><>< ======================================================================= ><>< */
/* ><><                        Sadiq, Hassan - SDQHAS001                        ><>< */
/* ><><                        PhD Statistics - STA6001W                        ><>< */
/* ><><                          ~~~~~~~~~~~~~~~~~~~~~                          ><>< */
/* ><><       ASSESSING THE ROBUSTNESS OF PHYLOGENETIC AND PHYLOGEOGRAPHIC      ><>< */
/* ><><            MODELS TO CHANGES IN RESIDUE FITNESSES OVER TIME:            ><>< */
/* ><><                           A Simulation Study.                           ><>< */
/* ><><                          ~~~~~~~~~~~~~~~~~~~~~                          ><>< */
/* ><><                               HyPhy Code:                               ><>< */
/* ><><          Illustrative Simulation for Literature Review Chapter          ><>< */
/* ><><                          ~~~~~~~~~~~~~~~~~~~~~                          ><>< */
/* ><><                             9th June, 2019.                             ><>< */
/* ><>< ======================================================================= ><>< */

/* ><>< ======================================================================= ><>< */
/* ><><                          COMPILATION RESOURCES                          ><>< */
/* ><><                          ~~~~~~~~~~~~~~~~~~~~~                          ><>< */
/*       Computer System (MacBook) :-                                                */
/*            Intel Core i7, CPU 2.90GHz, 16GB RAM                                   */
/*            MacBook Pro OS Mojave (Version 10.14.5)                                */
/*       HyPhy :-                                                                    */
/*            Version 2.3.1120180430 Beta(MPI) for Darwin on x86_64                  */
/* ><>< ======================================================================= ><>< */

k := 1e-06;
w := 1e-02;
eqAAfq = {61,1}["1/61"];
fprintf("hyphyTimes.txt", CLEAR_FILE);
siteSizes = {{5e+01, 5e+02, 5e+03, 5e+04, 5e+05}};
GY94 = {61, 61, {0,1,t*w}
                {0,2,k*t}     {0,3,t*w}     {0,4,t*w}   {0,8,k*t*w}    {0,12,t*w}
               {0,16,t*w}  {0,32,k*t*w}     {1,0,t*w}     {1,2,t*w}     {1,3,k*t}
                {1,5,t*w}   {1,9,k*t*w}    {1,13,t*w}    {1,17,t*w}  {1,33,k*t*w}
               {1,48,t*w}     {2,0,k*t}     {2,1,t*w}     {2,3,t*w}     {2,6,t*w}
             {2,10,k*t*w}    {2,14,t*w}    {2,18,t*w}  {2,34,k*t*w}     {3,0,t*w}
                {3,1,k*t}     {3,2,t*w}     {3,7,t*w}  {3,11,k*t*w}    {3,15,t*w}
               {3,19,t*w}  {3,35,k*t*w}    {3,49,t*w}     {4,0,t*w}       {4,5,t}
                {4,6,k*t}       {4,7,t}     {4,8,t*w}  {4,12,k*t*w}    {4,20,t*w}
             {4,36,k*t*w}    {4,50,t*w}     {5,1,t*w}       {5,4,t}       {5,6,t}
                {5,7,k*t}     {5,9,t*w}  {5,13,k*t*w}    {5,21,t*w}  {5,37,k*t*w}
               {5,51,t*w}     {6,2,t*w}     {6,4,k*t}       {6,5,t}       {6,7,t}
               {6,10,t*w}  {6,14,k*t*w}    {6,22,t*w}  {6,38,k*t*w}    {6,52,t*w}
                {7,3,t*w}       {7,4,t}     {7,5,k*t}       {7,6,t}    {7,11,t*w}
             {7,15,k*t*w}    {7,23,t*w}  {7,39,k*t*w}    {7,53,t*w}   {8,0,k*t*w}
                {8,4,t*w}     {8,9,t*w}    {8,10,k*t}    {8,11,t*w}    {8,12,t*w}
                 {8,24,t}  {8,40,k*t*w}   {9,1,k*t*w}     {9,5,t*w}     {9,8,t*w}
               {9,10,t*w}    {9,11,k*t}    {9,13,t*w}    {9,25,t*w}  {9,41,k*t*w}
               {9,54,t*w}  {10,2,k*t*w}    {10,6,t*w}    {10,8,k*t}    {10,9,t*w}
              {10,11,t*w}   {10,14,t*w}     {10,26,t} {10,42,k*t*w}   {10,55,t*w}
             {11,3,k*t*w}    {11,7,t*w}    {11,8,t*w}    {11,9,k*t}   {11,10,t*w}
              {11,15,t*w}   {11,27,t*w} {11,43,k*t*w}   {11,56,t*w}    {12,0,t*w}
             {12,4,k*t*w}    {12,8,t*w}     {12,13,t} {12,14,k*t*w}     {12,15,t}
              {12,28,t*w} {12,44,k*t*w}   {12,57,t*w}    {13,1,t*w}  {13,5,k*t*w}
               {13,9,t*w}     {13,12,t}   {13,14,t*w}   {13,15,k*t}   {13,29,t*w}
            {13,45,k*t*w}   {13,58,t*w}    {14,2,t*w}  {14,6,k*t*w}   {14,10,t*w}
            {14,12,k*t*w}   {14,13,t*w}   {14,15,t*w}   {14,30,t*w} {14,46,k*t*w}
              {14,59,t*w}    {15,3,t*w}  {15,7,k*t*w}   {15,11,t*w}     {15,12,t}
              {15,13,k*t}   {15,14,t*w}   {15,31,t*w} {15,47,k*t*w}   {15,60,t*w}
               {16,0,t*w}   {16,17,t*w}   {16,18,k*t}   {16,19,t*w}   {16,20,t*w}
            {16,24,k*t*w}   {16,28,t*w}   {16,32,t*w}    {17,1,t*w}   {17,16,t*w}
              {17,18,t*w}   {17,19,k*t}   {17,21,t*w} {17,25,k*t*w}   {17,29,t*w}
              {17,33,t*w} {17,48,k*t*w}    {18,2,t*w}   {18,16,k*t}   {18,17,t*w}
              {18,19,t*w}   {18,22,t*w} {18,26,k*t*w}   {18,30,t*w}   {18,34,t*w}
               {19,3,t*w}   {19,16,t*w}   {19,17,k*t}   {19,18,t*w}   {19,23,t*w}
            {19,27,k*t*w}   {19,31,t*w}   {19,35,t*w} {19,49,k*t*w}    {20,4,t*w}
              {20,16,t*w}     {20,21,t}   {20,22,k*t}     {20,23,t}   {20,24,t*w}
            {20,28,k*t*w}   {20,36,t*w} {20,50,k*t*w}    {21,5,t*w}   {21,17,t*w}
                {21,20,t}     {21,22,t}   {21,23,k*t}   {21,25,t*w} {21,29,k*t*w}
              {21,37,t*w} {21,51,k*t*w}    {22,6,t*w}   {22,18,t*w}   {22,20,k*t}
                {22,21,t}     {22,23,t}   {22,26,t*w} {22,30,k*t*w}   {22,38,t*w}
            {22,52,k*t*w}    {23,7,t*w}   {23,19,t*w}     {23,20,t}   {23,21,k*t}
                {23,22,t}   {23,27,t*w} {23,31,k*t*w}   {23,39,t*w} {23,53,k*t*w}
                 {24,8,t} {24,16,k*t*w}   {24,20,t*w}     {24,25,t}   {24,26,k*t}
                {24,27,t}   {24,28,t*w}   {24,40,t*w}    {25,9,t*w} {25,17,k*t*w}
              {25,21,t*w}     {25,24,t}     {25,26,t}   {25,27,k*t}   {25,29,t*w}
              {25,41,t*w} {25,54,k*t*w}     {26,10,t} {26,18,k*t*w}   {26,22,t*w}
              {26,24,k*t}     {26,25,t}     {26,27,t}   {26,30,t*w}   {26,42,t*w}
            {26,55,k*t*w}   {27,11,t*w} {27,19,k*t*w}   {27,23,t*w}     {27,24,t}
              {27,25,k*t}     {27,26,t}   {27,31,t*w}   {27,43,t*w} {27,56,k*t*w}
              {28,12,t*w}   {28,16,t*w} {28,20,k*t*w}   {28,24,t*w}     {28,29,t}
              {28,30,k*t}     {28,31,t}   {28,44,t*w}   {28,57,k*t}   {29,13,t*w}
              {29,17,t*w} {29,21,k*t*w}   {29,25,t*w}     {29,28,t}     {29,30,t}
              {29,31,k*t}   {29,45,t*w} {29,58,k*t*w}   {30,14,t*w}   {30,18,t*w}
            {30,22,k*t*w}   {30,26,t*w}   {30,28,k*t}     {30,29,t}     {30,31,t}
              {30,46,t*w}   {30,59,k*t}   {31,15,t*w}   {31,19,t*w} {31,23,k*t*w}
              {31,27,t*w}     {31,28,t}   {31,29,k*t}     {31,30,t}   {31,47,t*w}
            {31,60,k*t*w}  {32,0,k*t*w}   {32,16,t*w}   {32,33,t*w}   {32,34,k*t}
              {32,35,t*w}   {32,36,t*w} {32,40,k*t*w}   {32,44,t*w}  {33,1,k*t*w}
              {33,17,t*w}   {33,32,t*w}   {33,34,t*w}   {33,35,k*t}   {33,37,t*w}
            {33,41,k*t*w}   {33,45,t*w}   {33,48,t*w}  {34,2,k*t*w}   {34,18,t*w}
              {34,32,k*t}   {34,33,t*w}   {34,35,t*w}   {34,38,t*w} {34,42,k*t*w}
              {34,46,t*w}  {35,3,k*t*w}   {35,19,t*w}   {35,32,t*w}   {35,33,k*t}
              {35,34,t*w}   {35,39,t*w} {35,43,k*t*w}   {35,47,t*w}   {35,49,t*w}
             {36,4,k*t*w}   {36,20,t*w}   {36,32,t*w}     {36,37,t}   {36,38,k*t}
                {36,39,t}   {36,40,t*w} {36,44,k*t*w}   {36,50,t*w}  {37,5,k*t*w}
              {37,21,t*w}   {37,33,t*w}     {37,36,t}     {37,38,t}   {37,39,k*t}
              {37,41,t*w} {37,45,k*t*w}   {37,51,t*w}  {38,6,k*t*w}   {38,22,t*w}
              {38,34,t*w}   {38,36,k*t}     {38,37,t}     {38,39,t}   {38,42,t*w}
            {38,46,k*t*w}   {38,52,t*w}  {39,7,k*t*w}   {39,23,t*w}   {39,35,t*w}
                {39,36,t}   {39,37,k*t}     {39,38,t}   {39,43,t*w} {39,47,k*t*w}
              {39,53,t*w}  {40,8,k*t*w}   {40,24,t*w} {40,32,k*t*w}   {40,36,t*w}
                {40,41,t}   {40,42,k*t}     {40,43,t}   {40,44,t*w}  {41,9,k*t*w}
              {41,25,t*w} {41,33,k*t*w}   {41,37,t*w}     {41,40,t}     {41,42,t}
              {41,43,k*t}   {41,45,t*w}   {41,54,t*w} {42,10,k*t*w}   {42,26,t*w}
            {42,34,k*t*w}   {42,38,t*w}   {42,40,k*t}     {42,41,t}     {42,43,t}
              {42,46,t*w}   {42,55,t*w} {43,11,k*t*w}   {43,27,t*w} {43,35,k*t*w}
              {43,39,t*w}     {43,40,t}   {43,41,k*t}     {43,42,t}   {43,47,t*w}
              {43,56,t*w} {44,12,k*t*w}   {44,28,t*w}   {44,32,t*w} {44,36,k*t*w}
              {44,40,t*w}     {44,45,t}   {44,46,k*t}     {44,47,t}   {44,57,t*w}
            {45,13,k*t*w}   {45,29,t*w}   {45,33,t*w} {45,37,k*t*w}   {45,41,t*w}
                {45,44,t}     {45,46,t}   {45,47,k*t}   {45,58,t*w} {46,14,k*t*w}
              {46,30,t*w}   {46,34,t*w} {46,38,k*t*w}   {46,42,t*w}   {46,44,k*t}
                {46,45,t}     {46,47,t}   {46,59,t*w} {47,15,k*t*w}   {47,31,t*w}
              {47,35,t*w} {47,39,k*t*w}   {47,43,t*w}     {47,44,t}   {47,45,k*t}
                {47,46,t}   {47,60,t*w}    {48,1,t*w} {48,17,k*t*w}   {48,33,t*w}
              {48,49,k*t}   {48,51,t*w} {48,54,k*t*w}   {48,58,t*w}    {49,3,t*w}
            {49,19,k*t*w}   {49,35,t*w}   {49,48,k*t}   {49,53,t*w} {49,56,k*t*w}
              {49,60,t*w}    {50,4,t*w} {50,20,k*t*w}   {50,36,t*w}     {50,51,t}
              {50,52,k*t}     {50,53,t} {50,57,k*t*w}    {51,5,t*w} {51,21,k*t*w}
              {51,37,t*w}   {51,48,t*w}     {51,50,t}     {51,52,t}   {51,53,k*t}
              {51,54,t*w} {51,58,k*t*w}    {52,6,t*w} {52,22,k*t*w}   {52,38,t*w}
              {52,50,k*t}     {52,51,t}     {52,53,t}   {52,55,t*w} {52,59,k*t*w}
               {53,7,t*w} {53,23,k*t*w}   {53,39,t*w}   {53,49,t*w}     {53,50,t}
              {53,51,k*t}     {53,52,t}   {53,56,t*w} {53,60,k*t*w}    {54,9,t*w}
            {54,25,k*t*w}   {54,41,t*w} {54,48,k*t*w}   {54,51,t*w}   {54,55,t*w}
              {54,56,k*t}   {54,58,t*w}   {55,10,t*w} {55,26,k*t*w}   {55,42,t*w}
              {55,52,t*w}   {55,54,t*w}   {55,56,t*w}   {55,59,t*w}   {56,11,t*w}
            {56,27,k*t*w}   {56,43,t*w} {56,49,k*t*w}   {56,53,t*w}   {56,54,k*t}
              {56,55,t*w}   {56,60,t*w}   {57,12,t*w}   {57,28,k*t}   {57,44,t*w}
            {57,50,k*t*w}   {57,58,t*w}   {57,59,k*t}   {57,60,t*w}   {58,13,t*w}
            {58,29,k*t*w}   {58,45,t*w}   {58,48,t*w} {58,51,k*t*w}   {58,54,t*w}
              {58,57,t*w}   {58,59,t*w}   {58,60,k*t}   {59,14,t*w}   {59,30,k*t}
              {59,46,t*w} {59,52,k*t*w}   {59,55,t*w}   {59,57,k*t}   {59,58,t*w}
              {59,60,t*w}   {60,15,t*w} {60,31,k*t*w}   {60,47,t*w}   {60,49,t*w}
            {60,53,k*t*w}   {60,56,t*w}   {60,57,t*w}   {60,58,k*t}   {60,59,t*w}};
treeText = "(((((((((((Species0001:0.05,Species0002:0.05)iNode0:0.05,(Specie" +
  "s0003:0.05,Species0004:0.05)iNode1:0.05)iNode512:0.05,((Species0005:0.05,Speci" +
  "es0006:0.05)iNode2:0.05,(Species0007:0.05,Species0008:0.05)iNode3:0.05)iNode51" +
  "3:0.05)iNode768:0.05,(((Species0009:0.05,Species0010:0.05)iNode4:0.05,(Species" +
  "0011:0.05,Species0012:0.05)iNode5:0.05)iNode514:0.05,((Species0013:0.05,Specie" +
  "s0014:0.05)iNode6:0.05,(Species0015:0.05,Species0016:0.05)iNode7:0.05)iNode515" +
  ":0.05)iNode769:0.05)iNode896:0.05,((((Species0017:0.05,Species0018:0.05)iNode8" +
  ":0.05,(Species0019:0.05,Species0020:0.05)iNode9:0.05)iNode516:0.05,((Species00" +
  "21:0.05,Species0022:0.05)iNode10:0.05,(Species0023:0.05,Species0024:0.05)iNode" +
  "11:0.05)iNode517:0.05)iNode770:0.05,(((Species0025:0.05,Species0026:0.05)iNode" +
  "12:0.05,(Species0027:0.05,Species0028:0.05)iNode13:0.05)iNode518:0.05,((Specie" +
  "s0029:0.05,Species0030:0.05)iNode14:0.05,(Species0031:0.05,Species0032:0.05)iN" +
  "ode15:0.05)iNode519:0.05)iNode771:0.05)iNode897:0.05)iNode960:0.05,(((((Specie" +
  "s0033:0.05,Species0034:0.05)iNode16:0.05,(Species0035:0.05,Species0036:0.05)iN" +
  "ode17:0.05)iNode520:0.05,((Species0037:0.05,Species0038:0.05)iNode18:0.05,(Spe" +
  "cies0039:0.05,Species0040:0.05)iNode19:0.05)iNode521:0.05)iNode772:0.05,(((Spe" +
  "cies0041:0.05,Species0042:0.05)iNode20:0.05,(Species0043:0.05,Species0044:0.05" +
  ")iNode21:0.05)iNode522:0.05,((Species0045:0.05,Species0046:0.05)iNode22:0.05,(" +
  "Species0047:0.05,Species0048:0.05)iNode23:0.05)iNode523:0.05)iNode773:0.05)iNo" +
  "de898:0.05,((((Species0049:0.05,Species0050:0.05)iNode24:0.05,(Species0051:0.0" +
  "5,Species0052:0.05)iNode25:0.05)iNode524:0.05,((Species0053:0.05,Species0054:0" +
  ".05)iNode26:0.05,(Species0055:0.05,Species0056:0.05)iNode27:0.05)iNode525:0.05" +
  ")iNode774:0.05,(((Species0057:0.05,Species0058:0.05)iNode28:0.05,(Species0059:" +
  "0.05,Species0060:0.05)iNode29:0.05)iNode526:0.05,((Species0061:0.05,Species006" +
  "2:0.05)iNode30:0.05,(Species0063:0.05,Species0064:0.05)iNode31:0.05)iNode527:0" +
  ".05)iNode775:0.05)iNode899:0.05)iNode961:0.05)iNode992:0.05,((((((Species0065:" +
  "0.05,Species0066:0.05)iNode32:0.05,(Species0067:0.05,Species0068:0.05)iNode33:" +
  "0.05)iNode528:0.05,((Species0069:0.05,Species0070:0.05)iNode34:0.05,(Species00" +
  "71:0.05,Species0072:0.05)iNode35:0.05)iNode529:0.05)iNode776:0.05,(((Species00" +
  "73:0.05,Species0074:0.05)iNode36:0.05,(Species0075:0.05,Species0076:0.05)iNode" +
  "37:0.05)iNode530:0.05,((Species0077:0.05,Species0078:0.05)iNode38:0.05,(Specie" +
  "s0079:0.05,Species0080:0.05)iNode39:0.05)iNode531:0.05)iNode777:0.05)iNode900:" +
  "0.05,((((Species0081:0.05,Species0082:0.05)iNode40:0.05,(Species0083:0.05,Spec" +
  "ies0084:0.05)iNode41:0.05)iNode532:0.05,((Species0085:0.05,Species0086:0.05)iN" +
  "ode42:0.05,(Species0087:0.05,Species0088:0.05)iNode43:0.05)iNode533:0.05)iNode" +
  "778:0.05,(((Species0089:0.05,Species0090:0.05)iNode44:0.05,(Species0091:0.05,S" +
  "pecies0092:0.05)iNode45:0.05)iNode534:0.05,((Species0093:0.05,Species0094:0.05" +
  ")iNode46:0.05,(Species0095:0.05,Species0096:0.05)iNode47:0.05)iNode535:0.05)iN" +
  "ode779:0.05)iNode901:0.05)iNode962:0.05,(((((Species0097:0.05,Species0098:0.05" +
  ")iNode48:0.05,(Species0099:0.05,Species0100:0.05)iNode49:0.05)iNode536:0.05,((" +
  "Species0101:0.05,Species0102:0.05)iNode50:0.05,(Species0103:0.05,Species0104:0" +
  ".05)iNode51:0.05)iNode537:0.05)iNode780:0.05,(((Species0105:0.05,Species0106:0" +
  ".05)iNode52:0.05,(Species0107:0.05,Species0108:0.05)iNode53:0.05)iNode538:0.05" +
  ",((Species0109:0.05,Species0110:0.05)iNode54:0.05,(Species0111:0.05,Species011" +
  "2:0.05)iNode55:0.05)iNode539:0.05)iNode781:0.05)iNode902:0.05,((((Species0113:" +
  "0.05,Species0114:0.05)iNode56:0.05,(Species0115:0.05,Species0116:0.05)iNode57:" +
  "0.05)iNode540:0.05,((Species0117:0.05,Species0118:0.05)iNode58:0.05,(Species01" +
  "19:0.05,Species0120:0.05)iNode59:0.05)iNode541:0.05)iNode782:0.05,(((Species01" +
  "21:0.05,Species0122:0.05)iNode60:0.05,(Species0123:0.05,Species0124:0.05)iNode" +
  "61:0.05)iNode542:0.05,((Species0125:0.05,Species0126:0.05)iNode62:0.05,(Specie" +
  "s0127:0.05,Species0128:0.05)iNode63:0.05)iNode543:0.05)iNode783:0.05)iNode903:" +
  "0.05)iNode963:0.05)iNode993:0.05)iNode1008:0.05,(((((((Species0129:0.05,Specie" +
  "s0130:0.05)iNode64:0.05,(Species0131:0.05,Species0132:0.05)iNode65:0.05)iNode5" +
  "44:0.05,((Species0133:0.05,Species0134:0.05)iNode66:0.05,(Species0135:0.05,Spe" +
  "cies0136:0.05)iNode67:0.05)iNode545:0.05)iNode784:0.05,(((Species0137:0.05,Spe" +
  "cies0138:0.05)iNode68:0.05,(Species0139:0.05,Species0140:0.05)iNode69:0.05)iNo" +
  "de546:0.05,((Species0141:0.05,Species0142:0.05)iNode70:0.05,(Species0143:0.05," +
  "Species0144:0.05)iNode71:0.05)iNode547:0.05)iNode785:0.05)iNode904:0.05,((((Sp" +
  "ecies0145:0.05,Species0146:0.05)iNode72:0.05,(Species0147:0.05,Species0148:0.0" +
  "5)iNode73:0.05)iNode548:0.05,((Species0149:0.05,Species0150:0.05)iNode74:0.05," +
  "(Species0151:0.05,Species0152:0.05)iNode75:0.05)iNode549:0.05)iNode786:0.05,((" +
  "(Species0153:0.05,Species0154:0.05)iNode76:0.05,(Species0155:0.05,Species0156:" +
  "0.05)iNode77:0.05)iNode550:0.05,((Species0157:0.05,Species0158:0.05)iNode78:0." +
  "05,(Species0159:0.05,Species0160:0.05)iNode79:0.05)iNode551:0.05)iNode787:0.05" +
  ")iNode905:0.05)iNode964:0.05,(((((Species0161:0.05,Species0162:0.05)iNode80:0." +
  "05,(Species0163:0.05,Species0164:0.05)iNode81:0.05)iNode552:0.05,((Species0165" +
  ":0.05,Species0166:0.05)iNode82:0.05,(Species0167:0.05,Species0168:0.05)iNode83" +
  ":0.05)iNode553:0.05)iNode788:0.05,(((Species0169:0.05,Species0170:0.05)iNode84" +
  ":0.05,(Species0171:0.05,Species0172:0.05)iNode85:0.05)iNode554:0.05,((Species0" +
  "173:0.05,Species0174:0.05)iNode86:0.05,(Species0175:0.05,Species0176:0.05)iNod" +
  "e87:0.05)iNode555:0.05)iNode789:0.05)iNode906:0.05,((((Species0177:0.05,Specie" +
  "s0178:0.05)iNode88:0.05,(Species0179:0.05,Species0180:0.05)iNode89:0.05)iNode5" +
  "56:0.05,((Species0181:0.05,Species0182:0.05)iNode90:0.05,(Species0183:0.05,Spe" +
  "cies0184:0.05)iNode91:0.05)iNode557:0.05)iNode790:0.05,(((Species0185:0.05,Spe" +
  "cies0186:0.05)iNode92:0.05,(Species0187:0.05,Species0188:0.05)iNode93:0.05)iNo" +
  "de558:0.05,((Species0189:0.05,Species0190:0.05)iNode94:0.05,(Species0191:0.05," +
  "Species0192:0.05)iNode95:0.05)iNode559:0.05)iNode791:0.05)iNode907:0.05)iNode9" +
  "65:0.05)iNode994:0.05,((((((Species0193:0.05,Species0194:0.05)iNode96:0.05,(Sp" +
  "ecies0195:0.05,Species0196:0.05)iNode97:0.05)iNode560:0.05,((Species0197:0.05," +
  "Species0198:0.05)iNode98:0.05,(Species0199:0.05,Species0200:0.05)iNode99:0.05)" +
  "iNode561:0.05)iNode792:0.05,(((Species0201:0.05,Species0202:0.05)iNode100:0.05" +
  ",(Species0203:0.05,Species0204:0.05)iNode101:0.05)iNode562:0.05,((Species0205:" +
  "0.05,Species0206:0.05)iNode102:0.05,(Species0207:0.05,Species0208:0.05)iNode10" +
  "3:0.05)iNode563:0.05)iNode793:0.05)iNode908:0.05,((((Species0209:0.05,Species0" +
  "210:0.05)iNode104:0.05,(Species0211:0.05,Species0212:0.05)iNode105:0.05)iNode5" +
  "64:0.05,((Species0213:0.05,Species0214:0.05)iNode106:0.05,(Species0215:0.05,Sp" +
  "ecies0216:0.05)iNode107:0.05)iNode565:0.05)iNode794:0.05,(((Species0217:0.05,S" +
  "pecies0218:0.05)iNode108:0.05,(Species0219:0.05,Species0220:0.05)iNode109:0.05" +
  ")iNode566:0.05,((Species0221:0.05,Species0222:0.05)iNode110:0.05,(Species0223:" +
  "0.05,Species0224:0.05)iNode111:0.05)iNode567:0.05)iNode795:0.05)iNode909:0.05)" +
  "iNode966:0.05,(((((Species0225:0.05,Species0226:0.05)iNode112:0.05,(Species022" +
  "7:0.05,Species0228:0.05)iNode113:0.05)iNode568:0.05,((Species0229:0.05,Species" +
  "0230:0.05)iNode114:0.05,(Species0231:0.05,Species0232:0.05)iNode115:0.05)iNode" +
  "569:0.05)iNode796:0.05,(((Species0233:0.05,Species0234:0.05)iNode116:0.05,(Spe" +
  "cies0235:0.05,Species0236:0.05)iNode117:0.05)iNode570:0.05,((Species0237:0.05," +
  "Species0238:0.05)iNode118:0.05,(Species0239:0.05,Species0240:0.05)iNode119:0.0" +
  "5)iNode571:0.05)iNode797:0.05)iNode910:0.05,((((Species0241:0.05,Species0242:0" +
  ".05)iNode120:0.05,(Species0243:0.05,Species0244:0.05)iNode121:0.05)iNode572:0." +
  "05,((Species0245:0.05,Species0246:0.05)iNode122:0.05,(Species0247:0.05,Species" +
  "0248:0.05)iNode123:0.05)iNode573:0.05)iNode798:0.05,(((Species0249:0.05,Specie" +
  "s0250:0.05)iNode124:0.05,(Species0251:0.05,Species0252:0.05)iNode125:0.05)iNod" +
  "e574:0.05,((Species0253:0.05,Species0254:0.05)iNode126:0.05,(Species0255:0.05," +
  "Species0256:0.05)iNode127:0.05)iNode575:0.05)iNode799:0.05)iNode911:0.05)iNode" +
  "967:0.05)iNode995:0.05)iNode1009:0.05)iNode1016:0.05,((((((((Species0257:0.05," +
  "Species0258:0.05)iNode128:0.05,(Species0259:0.05,Species0260:0.05)iNode129:0.0" +
  "5)iNode576:0.05,((Species0261:0.05,Species0262:0.05)iNode130:0.05,(Species0263" +
  ":0.05,Species0264:0.05)iNode131:0.05)iNode577:0.05)iNode800:0.05,(((Species026" +
  "5:0.05,Species0266:0.05)iNode132:0.05,(Species0267:0.05,Species0268:0.05)iNode" +
  "133:0.05)iNode578:0.05,((Species0269:0.05,Species0270:0.05)iNode134:0.05,(Spec" +
  "ies0271:0.05,Species0272:0.05)iNode135:0.05)iNode579:0.05)iNode801:0.05)iNode9" +
  "12:0.05,((((Species0273:0.05,Species0274:0.05)iNode136:0.05,(Species0275:0.05," +
  "Species0276:0.05)iNode137:0.05)iNode580:0.05,((Species0277:0.05,Species0278:0." +
  "05)iNode138:0.05,(Species0279:0.05,Species0280:0.05)iNode139:0.05)iNode581:0.0" +
  "5)iNode802:0.05,(((Species0281:0.05,Species0282:0.05)iNode140:0.05,(Species028" +
  "3:0.05,Species0284:0.05)iNode141:0.05)iNode582:0.05,((Species0285:0.05,Species" +
  "0286:0.05)iNode142:0.05,(Species0287:0.05,Species0288:0.05)iNode143:0.05)iNode" +
  "583:0.05)iNode803:0.05)iNode913:0.05)iNode968:0.05,(((((Species0289:0.05,Speci" +
  "es0290:0.05)iNode144:0.05,(Species0291:0.05,Species0292:0.05)iNode145:0.05)iNo" +
  "de584:0.05,((Species0293:0.05,Species0294:0.05)iNode146:0.05,(Species0295:0.05" +
  ",Species0296:0.05)iNode147:0.05)iNode585:0.05)iNode804:0.05,(((Species0297:0.0" +
  "5,Species0298:0.05)iNode148:0.05,(Species0299:0.05,Species0300:0.05)iNode149:0" +
  ".05)iNode586:0.05,((Species0301:0.05,Species0302:0.05)iNode150:0.05,(Species03" +
  "03:0.05,Species0304:0.05)iNode151:0.05)iNode587:0.05)iNode805:0.05)iNode914:0." +
  "05,((((Species0305:0.05,Species0306:0.05)iNode152:0.05,(Species0307:0.05,Speci" +
  "es0308:0.05)iNode153:0.05)iNode588:0.05,((Species0309:0.05,Species0310:0.05)iN" +
  "ode154:0.05,(Species0311:0.05,Species0312:0.05)iNode155:0.05)iNode589:0.05)iNo" +
  "de806:0.05,(((Species0313:0.05,Species0314:0.05)iNode156:0.05,(Species0315:0.0" +
  "5,Species0316:0.05)iNode157:0.05)iNode590:0.05,((Species0317:0.05,Species0318:" +
  "0.05)iNode158:0.05,(Species0319:0.05,Species0320:0.05)iNode159:0.05)iNode591:0" +
  ".05)iNode807:0.05)iNode915:0.05)iNode969:0.05)iNode996:0.05,((((((Species0321:" +
  "0.05,Species0322:0.05)iNode160:0.05,(Species0323:0.05,Species0324:0.05)iNode16" +
  "1:0.05)iNode592:0.05,((Species0325:0.05,Species0326:0.05)iNode162:0.05,(Specie" +
  "s0327:0.05,Species0328:0.05)iNode163:0.05)iNode593:0.05)iNode808:0.05,(((Speci" +
  "es0329:0.05,Species0330:0.05)iNode164:0.05,(Species0331:0.05,Species0332:0.05)" +
  "iNode165:0.05)iNode594:0.05,((Species0333:0.05,Species0334:0.05)iNode166:0.05," +
  "(Species0335:0.05,Species0336:0.05)iNode167:0.05)iNode595:0.05)iNode809:0.05)i" +
  "Node916:0.05,((((Species0337:0.05,Species0338:0.05)iNode168:0.05,(Species0339:" +
  "0.05,Species0340:0.05)iNode169:0.05)iNode596:0.05,((Species0341:0.05,Species03" +
  "42:0.05)iNode170:0.05,(Species0343:0.05,Species0344:0.05)iNode171:0.05)iNode59" +
  "7:0.05)iNode810:0.05,(((Species0345:0.05,Species0346:0.05)iNode172:0.05,(Speci" +
  "es0347:0.05,Species0348:0.05)iNode173:0.05)iNode598:0.05,((Species0349:0.05,Sp" +
  "ecies0350:0.05)iNode174:0.05,(Species0351:0.05,Species0352:0.05)iNode175:0.05)" +
  "iNode599:0.05)iNode811:0.05)iNode917:0.05)iNode970:0.05,(((((Species0353:0.05," +
  "Species0354:0.05)iNode176:0.05,(Species0355:0.05,Species0356:0.05)iNode177:0.0" +
  "5)iNode600:0.05,((Species0357:0.05,Species0358:0.05)iNode178:0.05,(Species0359" +
  ":0.05,Species0360:0.05)iNode179:0.05)iNode601:0.05)iNode812:0.05,(((Species036" +
  "1:0.05,Species0362:0.05)iNode180:0.05,(Species0363:0.05,Species0364:0.05)iNode" +
  "181:0.05)iNode602:0.05,((Species0365:0.05,Species0366:0.05)iNode182:0.05,(Spec" +
  "ies0367:0.05,Species0368:0.05)iNode183:0.05)iNode603:0.05)iNode813:0.05)iNode9" +
  "18:0.05,((((Species0369:0.05,Species0370:0.05)iNode184:0.05,(Species0371:0.05," +
  "Species0372:0.05)iNode185:0.05)iNode604:0.05,((Species0373:0.05,Species0374:0." +
  "05)iNode186:0.05,(Species0375:0.05,Species0376:0.05)iNode187:0.05)iNode605:0.0" +
  "5)iNode814:0.05,(((Species0377:0.05,Species0378:0.05)iNode188:0.05,(Species037" +
  "9:0.05,Species0380:0.05)iNode189:0.05)iNode606:0.05,((Species0381:0.05,Species" +
  "0382:0.05)iNode190:0.05,(Species0383:0.05,Species0384:0.05)iNode191:0.05)iNode" +
  "607:0.05)iNode815:0.05)iNode919:0.05)iNode971:0.05)iNode997:0.05)iNode1010:0.0" +
  "5,(((((((Species0385:0.05,Species0386:0.05)iNode192:0.05,(Species0387:0.05,Spe" +
  "cies0388:0.05)iNode193:0.05)iNode608:0.05,((Species0389:0.05,Species0390:0.05)" +
  "iNode194:0.05,(Species0391:0.05,Species0392:0.05)iNode195:0.05)iNode609:0.05)i" +
  "Node816:0.05,(((Species0393:0.05,Species0394:0.05)iNode196:0.05,(Species0395:0" +
  ".05,Species0396:0.05)iNode197:0.05)iNode610:0.05,((Species0397:0.05,Species039" +
  "8:0.05)iNode198:0.05,(Species0399:0.05,Species0400:0.05)iNode199:0.05)iNode611" +
  ":0.05)iNode817:0.05)iNode920:0.05,((((Species0401:0.05,Species0402:0.05)iNode2" +
  "00:0.05,(Species0403:0.05,Species0404:0.05)iNode201:0.05)iNode612:0.05,((Speci" +
  "es0405:0.05,Species0406:0.05)iNode202:0.05,(Species0407:0.05,Species0408:0.05)" +
  "iNode203:0.05)iNode613:0.05)iNode818:0.05,(((Species0409:0.05,Species0410:0.05" +
  ")iNode204:0.05,(Species0411:0.05,Species0412:0.05)iNode205:0.05)iNode614:0.05," +
  "((Species0413:0.05,Species0414:0.05)iNode206:0.05,(Species0415:0.05,Species041" +
  "6:0.05)iNode207:0.05)iNode615:0.05)iNode819:0.05)iNode921:0.05)iNode972:0.05,(" +
  "((((Species0417:0.05,Species0418:0.05)iNode208:0.05,(Species0419:0.05,Species0" +
  "420:0.05)iNode209:0.05)iNode616:0.05,((Species0421:0.05,Species0422:0.05)iNode" +
  "210:0.05,(Species0423:0.05,Species0424:0.05)iNode211:0.05)iNode617:0.05)iNode8" +
  "20:0.05,(((Species0425:0.05,Species0426:0.05)iNode212:0.05,(Species0427:0.05,S" +
  "pecies0428:0.05)iNode213:0.05)iNode618:0.05,((Species0429:0.05,Species0430:0.0" +
  "5)iNode214:0.05,(Species0431:0.05,Species0432:0.05)iNode215:0.05)iNode619:0.05" +
  ")iNode821:0.05)iNode922:0.05,((((Species0433:0.05,Species0434:0.05)iNode216:0." +
  "05,(Species0435:0.05,Species0436:0.05)iNode217:0.05)iNode620:0.05,((Species043" +
  "7:0.05,Species0438:0.05)iNode218:0.05,(Species0439:0.05,Species0440:0.05)iNode" +
  "219:0.05)iNode621:0.05)iNode822:0.05,(((Species0441:0.05,Species0442:0.05)iNod" +
  "e220:0.05,(Species0443:0.05,Species0444:0.05)iNode221:0.05)iNode622:0.05,((Spe" +
  "cies0445:0.05,Species0446:0.05)iNode222:0.05,(Species0447:0.05,Species0448:0.0" +
  "5)iNode223:0.05)iNode623:0.05)iNode823:0.05)iNode923:0.05)iNode973:0.05)iNode9" +
  "98:0.05,((((((Species0449:0.05,Species0450:0.05)iNode224:0.05,(Species0451:0.0" +
  "5,Species0452:0.05)iNode225:0.05)iNode624:0.05,((Species0453:0.05,Species0454:" +
  "0.05)iNode226:0.05,(Species0455:0.05,Species0456:0.05)iNode227:0.05)iNode625:0" +
  ".05)iNode824:0.05,(((Species0457:0.05,Species0458:0.05)iNode228:0.05,(Species0" +
  "459:0.05,Species0460:0.05)iNode229:0.05)iNode626:0.05,((Species0461:0.05,Speci" +
  "es0462:0.05)iNode230:0.05,(Species0463:0.05,Species0464:0.05)iNode231:0.05)iNo" +
  "de627:0.05)iNode825:0.05)iNode924:0.05,((((Species0465:0.05,Species0466:0.05)i" +
  "Node232:0.05,(Species0467:0.05,Species0468:0.05)iNode233:0.05)iNode628:0.05,((" +
  "Species0469:0.05,Species0470:0.05)iNode234:0.05,(Species0471:0.05,Species0472:" +
  "0.05)iNode235:0.05)iNode629:0.05)iNode826:0.05,(((Species0473:0.05,Species0474" +
  ":0.05)iNode236:0.05,(Species0475:0.05,Species0476:0.05)iNode237:0.05)iNode630:" +
  "0.05,((Species0477:0.05,Species0478:0.05)iNode238:0.05,(Species0479:0.05,Speci" +
  "es0480:0.05)iNode239:0.05)iNode631:0.05)iNode827:0.05)iNode925:0.05)iNode974:0" +
  ".05,(((((Species0481:0.05,Species0482:0.05)iNode240:0.05,(Species0483:0.05,Spe" +
  "cies0484:0.05)iNode241:0.05)iNode632:0.05,((Species0485:0.05,Species0486:0.05)" +
  "iNode242:0.05,(Species0487:0.05,Species0488:0.05)iNode243:0.05)iNode633:0.05)i" +
  "Node828:0.05,(((Species0489:0.05,Species0490:0.05)iNode244:0.05,(Species0491:0" +
  ".05,Species0492:0.05)iNode245:0.05)iNode634:0.05,((Species0493:0.05,Species049" +
  "4:0.05)iNode246:0.05,(Species0495:0.05,Species0496:0.05)iNode247:0.05)iNode635" +
  ":0.05)iNode829:0.05)iNode926:0.05,((((Species0497:0.05,Species0498:0.05)iNode2" +
  "48:0.05,(Species0499:0.05,Species0500:0.05)iNode249:0.05)iNode636:0.05,((Speci" +
  "es0501:0.05,Species0502:0.05)iNode250:0.05,(Species0503:0.05,Species0504:0.05)" +
  "iNode251:0.05)iNode637:0.05)iNode830:0.05,(((Species0505:0.05,Species0506:0.05" +
  ")iNode252:0.05,(Species0507:0.05,Species0508:0.05)iNode253:0.05)iNode638:0.05," +
  "((Species0509:0.05,Species0510:0.05)iNode254:0.05,(Species0511:0.05,Species051" +
  "2:0.05)iNode255:0.05)iNode639:0.05)iNode831:0.05)iNode927:0.05)iNode975:0.05)i" +
  "Node999:0.05)iNode1011:0.05)iNode1017:0.05)iNode1020:0.05,(((((((((Species0513" +
  ":0.05,Species0514:0.05)iNode256:0.05,(Species0515:0.05,Species0516:0.05)iNode2" +
  "57:0.05)iNode640:0.05,((Species0517:0.05,Species0518:0.05)iNode258:0.05,(Speci" +
  "es0519:0.05,Species0520:0.05)iNode259:0.05)iNode641:0.05)iNode832:0.05,(((Spec" +
  "ies0521:0.05,Species0522:0.05)iNode260:0.05,(Species0523:0.05,Species0524:0.05" +
  ")iNode261:0.05)iNode642:0.05,((Species0525:0.05,Species0526:0.05)iNode262:0.05" +
  ",(Species0527:0.05,Species0528:0.05)iNode263:0.05)iNode643:0.05)iNode833:0.05)" +
  "iNode928:0.05,((((Species0529:0.05,Species0530:0.05)iNode264:0.05,(Species0531" +
  ":0.05,Species0532:0.05)iNode265:0.05)iNode644:0.05,((Species0533:0.05,Species0" +
  "534:0.05)iNode266:0.05,(Species0535:0.05,Species0536:0.05)iNode267:0.05)iNode6" +
  "45:0.05)iNode834:0.05,(((Species0537:0.05,Species0538:0.05)iNode268:0.05,(Spec" +
  "ies0539:0.05,Species0540:0.05)iNode269:0.05)iNode646:0.05,((Species0541:0.05,S" +
  "pecies0542:0.05)iNode270:0.05,(Species0543:0.05,Species0544:0.05)iNode271:0.05" +
  ")iNode647:0.05)iNode835:0.05)iNode929:0.05)iNode976:0.05,(((((Species0545:0.05" +
  ",Species0546:0.05)iNode272:0.05,(Species0547:0.05,Species0548:0.05)iNode273:0." +
  "05)iNode648:0.05,((Species0549:0.05,Species0550:0.05)iNode274:0.05,(Species055" +
  "1:0.05,Species0552:0.05)iNode275:0.05)iNode649:0.05)iNode836:0.05,(((Species05" +
  "53:0.05,Species0554:0.05)iNode276:0.05,(Species0555:0.05,Species0556:0.05)iNod" +
  "e277:0.05)iNode650:0.05,((Species0557:0.05,Species0558:0.05)iNode278:0.05,(Spe" +
  "cies0559:0.05,Species0560:0.05)iNode279:0.05)iNode651:0.05)iNode837:0.05)iNode" +
  "930:0.05,((((Species0561:0.05,Species0562:0.05)iNode280:0.05,(Species0563:0.05" +
  ",Species0564:0.05)iNode281:0.05)iNode652:0.05,((Species0565:0.05,Species0566:0" +
  ".05)iNode282:0.05,(Species0567:0.05,Species0568:0.05)iNode283:0.05)iNode653:0." +
  "05)iNode838:0.05,(((Species0569:0.05,Species0570:0.05)iNode284:0.05,(Species05" +
  "71:0.05,Species0572:0.05)iNode285:0.05)iNode654:0.05,((Species0573:0.05,Specie" +
  "s0574:0.05)iNode286:0.05,(Species0575:0.05,Species0576:0.05)iNode287:0.05)iNod" +
  "e655:0.05)iNode839:0.05)iNode931:0.05)iNode977:0.05)iNode1000:0.05,((((((Speci" +
  "es0577:0.05,Species0578:0.05)iNode288:0.05,(Species0579:0.05,Species0580:0.05)" +
  "iNode289:0.05)iNode656:0.05,((Species0581:0.05,Species0582:0.05)iNode290:0.05," +
  "(Species0583:0.05,Species0584:0.05)iNode291:0.05)iNode657:0.05)iNode840:0.05,(" +
  "((Species0585:0.05,Species0586:0.05)iNode292:0.05,(Species0587:0.05,Species058" +
  "8:0.05)iNode293:0.05)iNode658:0.05,((Species0589:0.05,Species0590:0.05)iNode29" +
  "4:0.05,(Species0591:0.05,Species0592:0.05)iNode295:0.05)iNode659:0.05)iNode841" +
  ":0.05)iNode932:0.05,((((Species0593:0.05,Species0594:0.05)iNode296:0.05,(Speci" +
  "es0595:0.05,Species0596:0.05)iNode297:0.05)iNode660:0.05,((Species0597:0.05,Sp" +
  "ecies0598:0.05)iNode298:0.05,(Species0599:0.05,Species0600:0.05)iNode299:0.05)" +
  "iNode661:0.05)iNode842:0.05,(((Species0601:0.05,Species0602:0.05)iNode300:0.05" +
  ",(Species0603:0.05,Species0604:0.05)iNode301:0.05)iNode662:0.05,((Species0605:" +
  "0.05,Species0606:0.05)iNode302:0.05,(Species0607:0.05,Species0608:0.05)iNode30" +
  "3:0.05)iNode663:0.05)iNode843:0.05)iNode933:0.05)iNode978:0.05,(((((Species060" +
  "9:0.05,Species0610:0.05)iNode304:0.05,(Species0611:0.05,Species0612:0.05)iNode" +
  "305:0.05)iNode664:0.05,((Species0613:0.05,Species0614:0.05)iNode306:0.05,(Spec" +
  "ies0615:0.05,Species0616:0.05)iNode307:0.05)iNode665:0.05)iNode844:0.05,(((Spe" +
  "cies0617:0.05,Species0618:0.05)iNode308:0.05,(Species0619:0.05,Species0620:0.0" +
  "5)iNode309:0.05)iNode666:0.05,((Species0621:0.05,Species0622:0.05)iNode310:0.0" +
  "5,(Species0623:0.05,Species0624:0.05)iNode311:0.05)iNode667:0.05)iNode845:0.05" +
  ")iNode934:0.05,((((Species0625:0.05,Species0626:0.05)iNode312:0.05,(Species062" +
  "7:0.05,Species0628:0.05)iNode313:0.05)iNode668:0.05,((Species0629:0.05,Species" +
  "0630:0.05)iNode314:0.05,(Species0631:0.05,Species0632:0.05)iNode315:0.05)iNode" +
  "669:0.05)iNode846:0.05,(((Species0633:0.05,Species0634:0.05)iNode316:0.05,(Spe" +
  "cies0635:0.05,Species0636:0.05)iNode317:0.05)iNode670:0.05,((Species0637:0.05," +
  "Species0638:0.05)iNode318:0.05,(Species0639:0.05,Species0640:0.05)iNode319:0.0" +
  "5)iNode671:0.05)iNode847:0.05)iNode935:0.05)iNode979:0.05)iNode1001:0.05)iNode" +
  "1012:0.05,(((((((Species0641:0.05,Species0642:0.05)iNode320:0.05,(Species0643:" +
  "0.05,Species0644:0.05)iNode321:0.05)iNode672:0.05,((Species0645:0.05,Species06" +
  "46:0.05)iNode322:0.05,(Species0647:0.05,Species0648:0.05)iNode323:0.05)iNode67" +
  "3:0.05)iNode848:0.05,(((Species0649:0.05,Species0650:0.05)iNode324:0.05,(Speci" +
  "es0651:0.05,Species0652:0.05)iNode325:0.05)iNode674:0.05,((Species0653:0.05,Sp" +
  "ecies0654:0.05)iNode326:0.05,(Species0655:0.05,Species0656:0.05)iNode327:0.05)" +
  "iNode675:0.05)iNode849:0.05)iNode936:0.05,((((Species0657:0.05,Species0658:0.0" +
  "5)iNode328:0.05,(Species0659:0.05,Species0660:0.05)iNode329:0.05)iNode676:0.05" +
  ",((Species0661:0.05,Species0662:0.05)iNode330:0.05,(Species0663:0.05,Species06" +
  "64:0.05)iNode331:0.05)iNode677:0.05)iNode850:0.05,(((Species0665:0.05,Species0" +
  "666:0.05)iNode332:0.05,(Species0667:0.05,Species0668:0.05)iNode333:0.05)iNode6" +
  "78:0.05,((Species0669:0.05,Species0670:0.05)iNode334:0.05,(Species0671:0.05,Sp" +
  "ecies0672:0.05)iNode335:0.05)iNode679:0.05)iNode851:0.05)iNode937:0.05)iNode98" +
  "0:0.05,(((((Species0673:0.05,Species0674:0.05)iNode336:0.05,(Species0675:0.05," +
  "Species0676:0.05)iNode337:0.05)iNode680:0.05,((Species0677:0.05,Species0678:0." +
  "05)iNode338:0.05,(Species0679:0.05,Species0680:0.05)iNode339:0.05)iNode681:0.0" +
  "5)iNode852:0.05,(((Species0681:0.05,Species0682:0.05)iNode340:0.05,(Species068" +
  "3:0.05,Species0684:0.05)iNode341:0.05)iNode682:0.05,((Species0685:0.05,Species" +
  "0686:0.05)iNode342:0.05,(Species0687:0.05,Species0688:0.05)iNode343:0.05)iNode" +
  "683:0.05)iNode853:0.05)iNode938:0.05,((((Species0689:0.05,Species0690:0.05)iNo" +
  "de344:0.05,(Species0691:0.05,Species0692:0.05)iNode345:0.05)iNode684:0.05,((Sp" +
  "ecies0693:0.05,Species0694:0.05)iNode346:0.05,(Species0695:0.05,Species0696:0." +
  "05)iNode347:0.05)iNode685:0.05)iNode854:0.05,(((Species0697:0.05,Species0698:0" +
  ".05)iNode348:0.05,(Species0699:0.05,Species0700:0.05)iNode349:0.05)iNode686:0." +
  "05,((Species0701:0.05,Species0702:0.05)iNode350:0.05,(Species0703:0.05,Species" +
  "0704:0.05)iNode351:0.05)iNode687:0.05)iNode855:0.05)iNode939:0.05)iNode981:0.0" +
  "5)iNode1002:0.05,((((((Species0705:0.05,Species0706:0.05)iNode352:0.05,(Specie" +
  "s0707:0.05,Species0708:0.05)iNode353:0.05)iNode688:0.05,((Species0709:0.05,Spe" +
  "cies0710:0.05)iNode354:0.05,(Species0711:0.05,Species0712:0.05)iNode355:0.05)i" +
  "Node689:0.05)iNode856:0.05,(((Species0713:0.05,Species0714:0.05)iNode356:0.05," +
  "(Species0715:0.05,Species0716:0.05)iNode357:0.05)iNode690:0.05,((Species0717:0" +
  ".05,Species0718:0.05)iNode358:0.05,(Species0719:0.05,Species0720:0.05)iNode359" +
  ":0.05)iNode691:0.05)iNode857:0.05)iNode940:0.05,((((Species0721:0.05,Species07" +
  "22:0.05)iNode360:0.05,(Species0723:0.05,Species0724:0.05)iNode361:0.05)iNode69" +
  "2:0.05,((Species0725:0.05,Species0726:0.05)iNode362:0.05,(Species0727:0.05,Spe" +
  "cies0728:0.05)iNode363:0.05)iNode693:0.05)iNode858:0.05,(((Species0729:0.05,Sp" +
  "ecies0730:0.05)iNode364:0.05,(Species0731:0.05,Species0732:0.05)iNode365:0.05)" +
  "iNode694:0.05,((Species0733:0.05,Species0734:0.05)iNode366:0.05,(Species0735:0" +
  ".05,Species0736:0.05)iNode367:0.05)iNode695:0.05)iNode859:0.05)iNode941:0.05)i" +
  "Node982:0.05,(((((Species0737:0.05,Species0738:0.05)iNode368:0.05,(Species0739" +
  ":0.05,Species0740:0.05)iNode369:0.05)iNode696:0.05,((Species0741:0.05,Species0" +
  "742:0.05)iNode370:0.05,(Species0743:0.05,Species0744:0.05)iNode371:0.05)iNode6" +
  "97:0.05)iNode860:0.05,(((Species0745:0.05,Species0746:0.05)iNode372:0.05,(Spec" +
  "ies0747:0.05,Species0748:0.05)iNode373:0.05)iNode698:0.05,((Species0749:0.05,S" +
  "pecies0750:0.05)iNode374:0.05,(Species0751:0.05,Species0752:0.05)iNode375:0.05" +
  ")iNode699:0.05)iNode861:0.05)iNode942:0.05,((((Species0753:0.05,Species0754:0." +
  "05)iNode376:0.05,(Species0755:0.05,Species0756:0.05)iNode377:0.05)iNode700:0.0" +
  "5,((Species0757:0.05,Species0758:0.05)iNode378:0.05,(Species0759:0.05,Species0" +
  "760:0.05)iNode379:0.05)iNode701:0.05)iNode862:0.05,(((Species0761:0.05,Species" +
  "0762:0.05)iNode380:0.05,(Species0763:0.05,Species0764:0.05)iNode381:0.05)iNode" +
  "702:0.05,((Species0765:0.05,Species0766:0.05)iNode382:0.05,(Species0767:0.05,S" +
  "pecies0768:0.05)iNode383:0.05)iNode703:0.05)iNode863:0.05)iNode943:0.05)iNode9" +
  "83:0.05)iNode1003:0.05)iNode1013:0.05)iNode1018:0.05,((((((((Species0769:0.05," +
  "Species0770:0.05)iNode384:0.05,(Species0771:0.05,Species0772:0.05)iNode385:0.0" +
  "5)iNode704:0.05,((Species0773:0.05,Species0774:0.05)iNode386:0.05,(Species0775" +
  ":0.05,Species0776:0.05)iNode387:0.05)iNode705:0.05)iNode864:0.05,(((Species077" +
  "7:0.05,Species0778:0.05)iNode388:0.05,(Species0779:0.05,Species0780:0.05)iNode" +
  "389:0.05)iNode706:0.05,((Species0781:0.05,Species0782:0.05)iNode390:0.05,(Spec" +
  "ies0783:0.05,Species0784:0.05)iNode391:0.05)iNode707:0.05)iNode865:0.05)iNode9" +
  "44:0.05,((((Species0785:0.05,Species0786:0.05)iNode392:0.05,(Species0787:0.05," +
  "Species0788:0.05)iNode393:0.05)iNode708:0.05,((Species0789:0.05,Species0790:0." +
  "05)iNode394:0.05,(Species0791:0.05,Species0792:0.05)iNode395:0.05)iNode709:0.0" +
  "5)iNode866:0.05,(((Species0793:0.05,Species0794:0.05)iNode396:0.05,(Species079" +
  "5:0.05,Species0796:0.05)iNode397:0.05)iNode710:0.05,((Species0797:0.05,Species" +
  "0798:0.05)iNode398:0.05,(Species0799:0.05,Species0800:0.05)iNode399:0.05)iNode" +
  "711:0.05)iNode867:0.05)iNode945:0.05)iNode984:0.05,(((((Species0801:0.05,Speci" +
  "es0802:0.05)iNode400:0.05,(Species0803:0.05,Species0804:0.05)iNode401:0.05)iNo" +
  "de712:0.05,((Species0805:0.05,Species0806:0.05)iNode402:0.05,(Species0807:0.05" +
  ",Species0808:0.05)iNode403:0.05)iNode713:0.05)iNode868:0.05,(((Species0809:0.0" +
  "5,Species0810:0.05)iNode404:0.05,(Species0811:0.05,Species0812:0.05)iNode405:0" +
  ".05)iNode714:0.05,((Species0813:0.05,Species0814:0.05)iNode406:0.05,(Species08" +
  "15:0.05,Species0816:0.05)iNode407:0.05)iNode715:0.05)iNode869:0.05)iNode946:0." +
  "05,((((Species0817:0.05,Species0818:0.05)iNode408:0.05,(Species0819:0.05,Speci" +
  "es0820:0.05)iNode409:0.05)iNode716:0.05,((Species0821:0.05,Species0822:0.05)iN" +
  "ode410:0.05,(Species0823:0.05,Species0824:0.05)iNode411:0.05)iNode717:0.05)iNo" +
  "de870:0.05,(((Species0825:0.05,Species0826:0.05)iNode412:0.05,(Species0827:0.0" +
  "5,Species0828:0.05)iNode413:0.05)iNode718:0.05,((Species0829:0.05,Species0830:" +
  "0.05)iNode414:0.05,(Species0831:0.05,Species0832:0.05)iNode415:0.05)iNode719:0" +
  ".05)iNode871:0.05)iNode947:0.05)iNode985:0.05)iNode1004:0.05,((((((Species0833" +
  ":0.05,Species0834:0.05)iNode416:0.05,(Species0835:0.05,Species0836:0.05)iNode4" +
  "17:0.05)iNode720:0.05,((Species0837:0.05,Species0838:0.05)iNode418:0.05,(Speci" +
  "es0839:0.05,Species0840:0.05)iNode419:0.05)iNode721:0.05)iNode872:0.05,(((Spec" +
  "ies0841:0.05,Species0842:0.05)iNode420:0.05,(Species0843:0.05,Species0844:0.05" +
  ")iNode421:0.05)iNode722:0.05,((Species0845:0.05,Species0846:0.05)iNode422:0.05" +
  ",(Species0847:0.05,Species0848:0.05)iNode423:0.05)iNode723:0.05)iNode873:0.05)" +
  "iNode948:0.05,((((Species0849:0.05,Species0850:0.05)iNode424:0.05,(Species0851" +
  ":0.05,Species0852:0.05)iNode425:0.05)iNode724:0.05,((Species0853:0.05,Species0" +
  "854:0.05)iNode426:0.05,(Species0855:0.05,Species0856:0.05)iNode427:0.05)iNode7" +
  "25:0.05)iNode874:0.05,(((Species0857:0.05,Species0858:0.05)iNode428:0.05,(Spec" +
  "ies0859:0.05,Species0860:0.05)iNode429:0.05)iNode726:0.05,((Species0861:0.05,S" +
  "pecies0862:0.05)iNode430:0.05,(Species0863:0.05,Species0864:0.05)iNode431:0.05" +
  ")iNode727:0.05)iNode875:0.05)iNode949:0.05)iNode986:0.05,(((((Species0865:0.05" +
  ",Species0866:0.05)iNode432:0.05,(Species0867:0.05,Species0868:0.05)iNode433:0." +
  "05)iNode728:0.05,((Species0869:0.05,Species0870:0.05)iNode434:0.05,(Species087" +
  "1:0.05,Species0872:0.05)iNode435:0.05)iNode729:0.05)iNode876:0.05,(((Species08" +
  "73:0.05,Species0874:0.05)iNode436:0.05,(Species0875:0.05,Species0876:0.05)iNod" +
  "e437:0.05)iNode730:0.05,((Species0877:0.05,Species0878:0.05)iNode438:0.05,(Spe" +
  "cies0879:0.05,Species0880:0.05)iNode439:0.05)iNode731:0.05)iNode877:0.05)iNode" +
  "950:0.05,((((Species0881:0.05,Species0882:0.05)iNode440:0.05,(Species0883:0.05" +
  ",Species0884:0.05)iNode441:0.05)iNode732:0.05,((Species0885:0.05,Species0886:0" +
  ".05)iNode442:0.05,(Species0887:0.05,Species0888:0.05)iNode443:0.05)iNode733:0." +
  "05)iNode878:0.05,(((Species0889:0.05,Species0890:0.05)iNode444:0.05,(Species08" +
  "91:0.05,Species0892:0.05)iNode445:0.05)iNode734:0.05,((Species0893:0.05,Specie" +
  "s0894:0.05)iNode446:0.05,(Species0895:0.05,Species0896:0.05)iNode447:0.05)iNod" +
  "e735:0.05)iNode879:0.05)iNode951:0.05)iNode987:0.05)iNode1005:0.05)iNode1014:0" +
  ".05,(((((((Species0897:0.05,Species0898:0.05)iNode448:0.05,(Species0899:0.05,S" +
  "pecies0900:0.05)iNode449:0.05)iNode736:0.05,((Species0901:0.05,Species0902:0.0" +
  "5)iNode450:0.05,(Species0903:0.05,Species0904:0.05)iNode451:0.05)iNode737:0.05" +
  ")iNode880:0.05,(((Species0905:0.05,Species0906:0.05)iNode452:0.05,(Species0907" +
  ":0.05,Species0908:0.05)iNode453:0.05)iNode738:0.05,((Species0909:0.05,Species0" +
  "910:0.05)iNode454:0.05,(Species0911:0.05,Species0912:0.05)iNode455:0.05)iNode7" +
  "39:0.05)iNode881:0.05)iNode952:0.05,((((Species0913:0.05,Species0914:0.05)iNod" +
  "e456:0.05,(Species0915:0.05,Species0916:0.05)iNode457:0.05)iNode740:0.05,((Spe" +
  "cies0917:0.05,Species0918:0.05)iNode458:0.05,(Species0919:0.05,Species0920:0.0" +
  "5)iNode459:0.05)iNode741:0.05)iNode882:0.05,(((Species0921:0.05,Species0922:0." +
  "05)iNode460:0.05,(Species0923:0.05,Species0924:0.05)iNode461:0.05)iNode742:0.0" +
  "5,((Species0925:0.05,Species0926:0.05)iNode462:0.05,(Species0927:0.05,Species0" +
  "928:0.05)iNode463:0.05)iNode743:0.05)iNode883:0.05)iNode953:0.05)iNode988:0.05" +
  ",(((((Species0929:0.05,Species0930:0.05)iNode464:0.05,(Species0931:0.05,Specie" +
  "s0932:0.05)iNode465:0.05)iNode744:0.05,((Species0933:0.05,Species0934:0.05)iNo" +
  "de466:0.05,(Species0935:0.05,Species0936:0.05)iNode467:0.05)iNode745:0.05)iNod" +
  "e884:0.05,(((Species0937:0.05,Species0938:0.05)iNode468:0.05,(Species0939:0.05" +
  ",Species0940:0.05)iNode469:0.05)iNode746:0.05,((Species0941:0.05,Species0942:0" +
  ".05)iNode470:0.05,(Species0943:0.05,Species0944:0.05)iNode471:0.05)iNode747:0." +
  "05)iNode885:0.05)iNode954:0.05,((((Species0945:0.05,Species0946:0.05)iNode472:" +
  "0.05,(Species0947:0.05,Species0948:0.05)iNode473:0.05)iNode748:0.05,((Species0" +
  "949:0.05,Species0950:0.05)iNode474:0.05,(Species0951:0.05,Species0952:0.05)iNo" +
  "de475:0.05)iNode749:0.05)iNode886:0.05,(((Species0953:0.05,Species0954:0.05)iN" +
  "ode476:0.05,(Species0955:0.05,Species0956:0.05)iNode477:0.05)iNode750:0.05,((S" +
  "pecies0957:0.05,Species0958:0.05)iNode478:0.05,(Species0959:0.05,Species0960:0" +
  ".05)iNode479:0.05)iNode751:0.05)iNode887:0.05)iNode955:0.05)iNode989:0.05)iNod" +
  "e1006:0.05,((((((Species0961:0.05,Species0962:0.05)iNode480:0.05,(Species0963:" +
  "0.05,Species0964:0.05)iNode481:0.05)iNode752:0.05,((Species0965:0.05,Species09" +
  "66:0.05)iNode482:0.05,(Species0967:0.05,Species0968:0.05)iNode483:0.05)iNode75" +
  "3:0.05)iNode888:0.05,(((Species0969:0.05,Species0970:0.05)iNode484:0.05,(Speci" +
  "es0971:0.05,Species0972:0.05)iNode485:0.05)iNode754:0.05,((Species0973:0.05,Sp" +
  "ecies0974:0.05)iNode486:0.05,(Species0975:0.05,Species0976:0.05)iNode487:0.05)" +
  "iNode755:0.05)iNode889:0.05)iNode956:0.05,((((Species0977:0.05,Species0978:0.0" +
  "5)iNode488:0.05,(Species0979:0.05,Species0980:0.05)iNode489:0.05)iNode756:0.05" +
  ",((Species0981:0.05,Species0982:0.05)iNode490:0.05,(Species0983:0.05,Species09" +
  "84:0.05)iNode491:0.05)iNode757:0.05)iNode890:0.05,(((Species0985:0.05,Species0" +
  "986:0.05)iNode492:0.05,(Species0987:0.05,Species0988:0.05)iNode493:0.05)iNode7" +
  "58:0.05,((Species0989:0.05,Species0990:0.05)iNode494:0.05,(Species0991:0.05,Sp" +
  "ecies0992:0.05)iNode495:0.05)iNode759:0.05)iNode891:0.05)iNode957:0.05)iNode99" +
  "0:0.05,(((((Species0993:0.05,Species0994:0.05)iNode496:0.05,(Species0995:0.05," +
  "Species0996:0.05)iNode497:0.05)iNode760:0.05,((Species0997:0.05,Species0998:0." +
  "05)iNode498:0.05,(Species0999:0.05,Species1000:0.05)iNode499:0.05)iNode761:0.0" +
  "5)iNode892:0.05,(((Species1001:0.05,Species1002:0.05)iNode500:0.05,(Species100" +
  "3:0.05,Species1004:0.05)iNode501:0.05)iNode762:0.05,((Species1005:0.05,Species" +
  "1006:0.05)iNode502:0.05,(Species1007:0.05,Species1008:0.05)iNode503:0.05)iNode" +
  "763:0.05)iNode893:0.05)iNode958:0.05,((((Species1009:0.05,Species1010:0.05)iNo" +
  "de504:0.05,(Species1011:0.05,Species1012:0.05)iNode505:0.05)iNode764:0.05,((Sp" +
  "ecies1013:0.05,Species1014:0.05)iNode506:0.05,(Species1015:0.05,Species1016:0." +
  "05)iNode507:0.05)iNode765:0.05)iNode894:0.05,(((Species1017:0.05,Species1018:0" +
  ".05)iNode508:0.05,(Species1019:0.05,Species1020:0.05)iNode509:0.05)iNode766:0." +
  "05,((Species1021:0.05,Species1022:0.05)iNode510:0.05,(Species1023:0.05,Species" +
  "1024:0.05)iNode511:0.05)iNode767:0.05)iNode895:0.05)iNode959:0.05)iNode991:0.0" +
  "5)iNode1007:0.05)iNode1015:0.05)iNode1019:0.05)iNode1021:0.05)iNode1022:0.05);";

for (i=0; i<Columns(siteSizes); i+=1){
  
  sites = siteSizes[i];
  
  for (j=1; j<=51; j+=1){
    startTime = Time(0);
    Model goldman = (GY94, eqAAfq, 1);
    Topology simulatedTopology = treeText;
    Tree simAtree = treeText;
    simXters = {{"A","C","G","T"}{"3","TAA,TAG,TGA","",""}};
    DataSet simulated = Simulate(simAtree, eqAAfq, simXters, sites, 0, "hyphySeqs.nex");
    cpuTime = Time(0) - startTime;
    fprintf("hyphyTimes.txt", cpuTime, "\t");
  }
  fprintf("hyphyTimes.txt", "\n");
  fprintf(stdout, "Completed simulations for size ", sites, "!.\n");
}

/* ><>< ======================================================================= ><>< */
/* ><><                             CODE ENDS HERE.                             ><>< */
/* ><>< ======================================================================= ><>< */