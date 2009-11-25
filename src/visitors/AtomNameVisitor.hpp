/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */

#ifndef VISITORS_ATOMNAMEVISITOR_HPP
#define VISITORS_ATOMNAMEVISITOR_HPP
#include "visitors/BaseVisitor.hpp"
#include "visitors/AtomData.hpp"
#include <string>
#include <map>
#include "visitors/AtomVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"

namespace OpenMD {


class AtomNameVisitor : public BaseVisitor {
    public:
      AtomNameVisitor(SimInfo* info);  
      AtomNameVisitor(SimInfo* info, const std::string& atomTypeFile);

      virtual void visit(Atom* atom) {visitAtom(atom);}
      virtual void visit(DirectionalAtom* datom) {visitAtom(static_cast<Atom*>(datom));}
      
      virtual void visit(RigidBody* rb);

      virtual const std::string toString();    
    private:
      void visitAtom(Atom* atom);
      void readAtomTypes(std::istream& is);
      std::string getBaseAtomTypeName(const std::string& atomTypeName);

      typedef std::map<std::string, std::string> MapType;
      MapType atomNames_;
      SimInfo* info_;
    };


const char defaultAtomTypeTable[] = { 
    "#Atom Type Base Atom Type  Atomic Number   Atomic Mass     Van derWaals Radius     Covalent Radius     Red Green   Blue \n "
    "X  X   0   0.0 1.0 0.0 255 20  147 \n "
    "H  H   1   1.008   1.2 0.32    250 235 215 \n "
    "He He  2   4.003   1.4 0.93    255 192 203 \n "
    "Li Li  3   6.941   1.82    1.23    178 34  34 \n "
    "Be Be  4   9.0122  1.3725  0.9 34  139 34 \n "
    "B  B   5   10.811  0.795   0.82    0   255 0 \n "
    "C  C   6   12.011  1.7 0.77    112 128 144 \n "
    "N  N   7   14.007  1.55    0.75    0   191 255 \n "
    "O  O   8   15.999  1.52    0.73    255 0   0 \n "
    "F  F   9   18.998  1.47    0.72    218 165 32 \n "
    "Ne Ne  10  20.18   1.54    0.71    255 105 180 \n "
    "Na Na  11  22.99   2.27    1.54    0   0   255 \n "
    "Na+     Na      11      22.99   2.27    1.54    0       0       255 \n "
    "Mg Mg  12  24.312  1.73    1.36    34  139 34 \n "
    "Al Al  13  26.982  1.7 1.18    190 190 190 \n "
    "Si Si  14  28.086  2.1 1.11    218 165 32 \n "
    "P  P   15  30.974  1.8 1.06    255 165 0 \n "
    "S  S   16  32.06   1.8 1.02    255 255 0 \n "
    "Cl Cl  17  35.453  1.75    0.99    0   255 0 \n "
    "Cl-      Cl      17      35.453  1.75    0.99    0       255     0 \n "
    "Ar Ar  18  39.948  1.88    0.98    255 192 203 \n "
    "K  K   19  39.098  2.75    2.03    255 20  147 \n "
    "Ca Ca  20  40.078  2.45    1.74    128 128 128 \n "
    "Sc Sc  21  44.956  1.37    1.44    190 190 190 \n "
    "Ti Ti  22  47.9    1.37    1.32    190 190 190 \n "
    "V  V   23  50.941  1.37    1.22    190 190 190 \n "
    "Cr Cr  24  51.996  1.37    1.18    190 190 190 \n "
    "Mn Mn  25  54.938  1.37    1.17    190 190 190 \n "
    "Fe Fe  26  55.847  1.456   1.17    255 165 0 \n "
    "Co Co  27  58.933  0.88    1.16    165 42  42 \n "
    "Ni Ni  28  58.71   0.69    1.15    165 42  42 \n "
    "Cu Cu  29  63.54   0.72    1.17    165 42  42 \n "
    "Zn Zn  30  65.37   0.74    1.25    165 42  42 \n "
    "Ga Ga  31  69.72   1.37    1.26    165 42  42 \n "
    "Ge Ge  32  72.59   1.95    1.22    85  107 47 \n "
    "As As  33  74.9216 1.85    1.2 253 245 230 \n "
    "Se Se  34  78.96   1.9 1.16    152 251 152 \n "
    "Br Br  35  79.904  1.85    1.14    165 42  42 \n "
    "Kr Kr  36  83.8    2.02    1.12    50  205 50 \n "
    "Rb Rb  37  85.47   1.58    2.16    165 42  42 \n "
    "Sr Sr  38  87.62   2.151   1.91    190 190 190 \n "
    "Y  Y   39  88.9059 1.801   1.62    190 190 190 \n "
    "Zr Zr  40  91.224  1.602   1.45    190 190 190 \n "
    "Nb Nb  41  92.9064 1.468   1.34    190 190 190 \n "
    "Mo Mo  42  95.94   1.526   1.3 255 127 80 \n "
    "Tc Tc  43  98.0    1.36    1.27    190 190 190 \n "
    "Ru Ru  44  101.07  1.339   1.25    190 190 190 \n "
    "Rh Rh  45  102.906 1.345   1.25    190 190 190 \n "
    "Pd Pd  46  106.42  1.376   1.28    190 190 190 \n "
    "Ag Ag  47  107.87  1.27    1.34    190 190 190 \n "
    "Cd Cd  48  112.41  1.424   1.48    255 140 0 \n "
    "In In  49  114.82  1.663   1.44    190 190 190 \n "
    "Sn Sn  50  118.71  2.1 1.41    190 190 190 \n "
    "Sb Sb  51  121.75  2.05    1.4 190 190 190 \n "
    "Te Te  52  127.6   2.06    1.36    190 190 190 \n "
    "I  I   53  129.905 1.98    1.33    160 32  240 \n "
    "Xe Xe  54  131.29  2.0 1.31    255 105 180 \n "
    "Cs Cs  55  132.905 1.84    2.35    165 42  42 \n "
    "Ba Ba  56  137.33  2.243   1.98    190 190 190 \n "
    "La La  57  138.906 1.877   1.69    190 190 190 \n "
    "Lu Lu  71  174.967 2.17    1.6 190 190 190 \n "
    "Hf Hf  72  178.49  1.58    1.44    190 190 190 \n "
    "Ta Ta  73  180.948 1.467   1.34    190 190 190 \n "
    "W  W   74  183.85  1.534   1.3 64  224 208 \n "
    "Re Re  75  186.207 1.375   1.28    190 190 190 \n "
    "Os Os  76  190.2   1.353   1.26    190 190 190 \n "
    "Ir Ir  77  192.22  1.357   1.27    190 190 190 \n "
    "Pt Pt  78  195.08  1.75    1.3 190 190 190 \n "
    "Au Au  79  196.967 1.66    1.34    255 215 0 \n "
    "Hg Hg  80  200.59  1.55    1.49    190 190 190 \n "
    "Tl Tl  81  204.383 1.96    1.48    190 190 190 \n "
    "Pb Pb  82  207.2   2.02    1.47    190 190 190 \n "
    "Bi Bi  83  208.98  2.15    1.46    255 181 197 \n "
    "LP X   0   0.0 1.0 0.0 255 20  147 \n "
    "PD X   0   0.0 0.0 0.0 0   0   0 \n "
    "DM X   0   0.0 0.0 0.0 0   0   0 \n "
    "1  H   1   1.008   1.2 0.32    250 235 215 \n "
    "H3 H   1   1.008   1.25    0.32    250 235 215 \n "
    "HALI   H   1   1.008   1.25    0.32    250 235 215 \n "
    "HARO   H   1   1.008   1.25    0.32    250 235 215 \n "
    "HC H   1   1.008   1.25    0.32    250 235 215 \n "
    "HHBN   H   1   1.008   1.25    0.32    250 235 215 \n "
    "HN H   1   1.008   1.25    0.32    250 235 215 \n "
    "HO H   1   1.008   1.25    0.32    250 235 215 \n "
    "HS H   1   1.008   1.25    0.32    250 235 215 \n "
    "HE He  2   4.003   1.4 0.93    255 192 203 \n "
    "3  Li  3   6.941   1.82    1.23    178 34  34 \n "
    "LI Li  3   6.941   1.82    1.23    178 34  34 \n "
    "5  B   5   10.811  0.795   0.82    0   255 0 \n "
    "6  C   6   12.011  1.7 0.77    112 128 144 \n "
    "C1 C   6   12.011  1.5 0.77    112 128 144 \n "
    "C2 C   6   12.011  1.5 0.77    112 128 144 \n "
    "C3 C   6   12.011  1.65    0.77    112 128 144 \n "
    "CA C   6   12.011  1.5 0.77    112 128 144 \n "
    "CAR    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CARO   C   6   12.011  1.5 0.77    112 128 144 \n "
    "CB C   6   12.011  1.5 0.77    112 128 144 \n "
    "CD C   6   12.011  1.5 0.77    112 128 144 \n "
    "CD1    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CD2    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CE1    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CE2    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CG C   6   12.011  1.5 0.77    112 128 144 \n "
    "CG1    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CG2    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CZ C   6   12.011  1.5 0.77    112 128 144 \n "
    "CH2AL  C   6   14.027  1.65    0.77    112 128 144 \n "
    "CH2OL  C   6   14.027  1.65    0.77    112 128 144 \n "
    "CH3AL  C   6   15.035  1.77    0.77    112 128 144 \n "
    "CHAL   C   6   13.019  1.67    0.77    112 128 144 \n "
    "CHAR   C   6   13.019  1.58    0.77    112 128 144 \n "
    "CHOL   C   6   13.019  1.67    0.77    112 128 144 \n "
    "CSP    C   6   12.011  1.5 0.77    112 128 144 \n "
    "CSP2   C   6   12.011  1.5 0.77    112 128 144 \n "
    "CSP3   C   6   12.011  1.65    0.77    112 128 144 \n "
    "CT C   6   15.035  1.77    0.77    112 128 144 \n "
    "7  N   7   14.007  1.55    0.75    0   191 255 \n "
    "NE N   7   14.007  1.35    0.75    0   191 255 \n "
    "ND2    N   7   14.007  1.35    0.75    0   191 255 \n "
    "NH1    N   7   14.007  1.35    0.75    0   191 255 \n "
    "NH2    N   7   14.007  1.35    0.75    0   191 255 \n "
    "N1 N   7   14.007  1.35    0.75    0   191 255 \n "
    "N2 N   7   14.007  1.35    0.75    0   191 255 \n "
    "N3 N   7   14.007  1.5 0.75    0   191 255 \n "
    "N3+    N   7   14.007  1.5 0.75    0   191 255 \n "
    "NAM    N   7   14.007  1.35    0.75    0   191 255 \n "
    "NAR    N   7   14.007  1.35    0.75    0   191 255 \n "
    "NARO   N   7   14.007  1.35    0.75    0   191 255 \n "
    "NPL3   N   7   14.007  1.35    0.75    0   191 255 \n "
    "NSP    N   7   14.007  1.35    0.75    0   191 255 \n "
    "NSP2   N   7   14.007  1.35    0.75    0   191 255 \n "
    "NSP3   N   7   14.007  1.5 0.75    0   191 255 \n "
    "NSP3+  N   7   14.007  1.5 0.75    0   191 255 \n "
    "8  O   8   15.999  1.52    0.73    255 0   0 \n "
    "O1      O       8       15.999  1.35    0.73    255     0       0 \n "
    "O2 O   8   15.999  1.35    0.73    255 0   0 \n "
    "O3 O   8   15.999  1.35    0.73    255 0   0 \n "
    "OH O   8   15.999  1.35    0.73    255 0   0 \n "
    "OD1    O   8   15.999  1.35    0.73    255 0   0 \n "
    "OE1    O   8   15.999  1.35    0.73    255 0   0 \n "
    "OE2    O   8   15.999  1.35    0.73    255 0   0 \n "
    "OG O   8   15.999  1.35    0.73    255 0   0 \n "
    "OG1    O   8   15.999  1.35    0.73    255 0   0 \n "
    "OXT    O   8   15.999  1.35    0.73    255 0   0 \n "
    "OSP2   O   8   15.999  1.35    0.73    255 0   0 \n "
    "OSP3   O   8   15.999  1.35    0.73    255 0   0 \n "
    "OSP3-  O   8   15.999  1.35    0.73    255 0   0 \n "
    "9  F   9   18.998  1.47    0.72    218 165 32 \n "
    "NE Ne  10  20.18   1.54    0.71    255 105 180 \n "
    "11 Na  11  22.99   2.27    1.54    0   0   255 \n "
    "NA Na  11  22.99   2.27    1.54    0   0   255 \n "
    "12 Mg  12  24.312  1.73    1.36    34  139 34 \n "
    "MG Mg  12  24.312  1.73    1.36    34  139 34 \n "
    "13 Al  13  26.982  1.7 1.18    190 190 190 \n "
    "AL Al  13  26.982  1.7 1.18    190 190 190 \n "
    "14 Si  14  28.086  2.1 1.11    218 165 32 \n "
    "SI Si  14  28.086  2.1 1.11    218 165 32 \n "
    "15 P   15  30.974  1.8 1.06    255 165 0 \n "
    "P(C)   P   15  30.974  1.75    1.06    255 165 0 \n "
    "P(N)   P   15  30.974  1.75    1.06    255 165 0 \n "
    "P(O)   P   15  30.974  1.75    1.06    255 165 0 \n "
    "P(S)   P   15  30.974  1.75    1.06    255 165 0 \n "
    "16 S   16  32.06   1.8 1.02    255 255 0 \n "
    "S2 S   16  32.06   1.85    1.02    255 255 0 \n "
    "S3 S   16  32.06   1.85    1.02    255 255 0 \n "
    "SG S   16  32.06   1.85    1.02    255 255 0 \n "
    "SDIV   S   16  32.06   1.85    1.02    255 255 0 \n "
    "STET   S   16  32.06   1.85    1.02    255 255 0 \n "
    "STRI   S   16  32.06   1.85    1.02    255 255 0 \n "
    "Sh S   16  32.06   1.85    1.02    255 255 0 \n "
    "17 Cl  17  35.453  1.75    0.99    0   255 0 \n "
    "CL Cl  17  35.453  1.75    0.99    0   255 0 \n "
    "AR Ar  18  39.948  1.88    0.98    255 192 203 \n "
    "HV Ar  18  399.48  1.91    0.0 255 192 203 \n "
    "LT Ar  18  39.948  1.91    0.0 255 192 203 \n "
    "DZ Ar  18  39.948  1.91    0.0 255 192 203 \n "
    "19 K   19  39.098  2.75    2.03    255 20  147 \n "
    "20 Ca  20  40.078  2.45    1.74    128 128 128 \n "
    "21 Sc  21  44.956  1.37    1.44    190 190 190 \n "
    "SC Sc  21  44.956  1.37    1.44    190 190 190 \n "
    "22 Ti  22  47.9    1.37    1.32    190 190 190 \n "
    "TI Ti  22  47.9    1.37    1.32    190 190 190 \n "
    "23 V   23  50.941  1.37    1.22    190 190 190 \n "
    "24 Cr  24  51.996  1.37    1.18    190 190 190 \n "
    "CR Cr  24  51.996  1.37    1.18    190 190 190 \n "
    "25 Mn  25  54.938  1.37    1.17    190 190 190 \n "
    "MN Mn  25  54.938  1.37    1.17    190 190 190 \n "
    "26 Fe  26  55.847  1.456   1.17    255 165 0 \n "
    "FE Fe  26  55.847  1.456   1.17    255 165 0 \n "
    "27 Co  27  58.933  0.88    1.16    165 42  42 \n "
    "CO Co  27  58.933  0.88    1.16    165 42  42 \n "
    "28 Ni  28  58.71   0.69    1.15    165 42  42 \n "
    "NI Ni  28  58.71   0.69    1.15    165 42  42 \n "
    "29 Cu  29  63.54   0.72    1.17    165 42  42 \n "
    "CU Cu  29  63.54   0.72    1.17    165 42  42 \n "
    "30 Zn  30  65.37   0.74    1.25    165 42  42 \n "
    "ZN Zn  30  65.37   0.74    1.25    165 42  42 \n "
    "31 Ga  31  69.72   1.37    1.26    165 42  42 \n "
    "GA Ga  31  69.72   1.37    1.26    165 42  42 \n "
    "GE Ge  32  72.59   1.95    1.22    85  107 47 \n "
    "35 Br  35  79.904  1.85    1.14    165 42  42 \n "
    "BR Br  35  79.904  1.85    1.14    165 42  42 \n "
    "37 Rb  37  85.47   1.58    2.16    165 42  42 \n "
    "RB Rb  37  85.47   1.58    2.16    165 42  42 \n "
    "MO Mo  42  95.94   1.526   1.3 255 127 80 \n "
    "47 Ag  47  107.87  1.27    1.34    190 190 190 \n "
    "AG Ag  47  107.87  1.27    1.34    190 190 190 \n "
    "53 I   53  129.905 1.98    1.33    160 32  240 \n "
    "55 Cs  55  132.905 1.84    2.35    165 42  42 \n "
    "CS Cs  55  132.905 1.84    2.35    165 42  42 \n "
    "79 Au  79  196.967 1.66    1.34    255 215 0 \n "
    "AU Au  79  196.967 1.66    1.34    255 215 0 \n "
    "XX X   0   0.0 1.0 0.0 255 20  147 \n "
    "+  X   0   0.0 1.0 1.38    255 20  147 \n "
    "++ X   0   0.0 1.0 0.9 255 20  147 \n "
    "-  X   0   0.0 1.0 2.0 255 20  147 \n "
    "-- X   0   0.0 1.0 2.0 255 20  147 \n "
    "DUMMY  X   0   0.0 1.0 0.0 255 20  147 \n "
    "RESERV X   0   0.0 1.0 0.0 255 20  147 \n "
    "BOGUS  X   0   0.0 1.0 0.0 255 20  147 \n "
    "DU X   0   0.0 1.0 0.0 255 20  147 \n "
    "Tv X   0   0.0 1.0 0.0 255 20  147 \n "
    "TV X   0   0.0 1.0 0.0 255 20  147 \n "
    "BQ X   0   0.0 1.0 0.0 255 20  147 \n "
    "GB X   0   760.09  1.0 0.0 255 20  147 \n "
    "GBDP   X   0   760.09  10.0    0.0 255 20  147 \n "
    "linear linear  0   760.09  1.0 0.0 255 20  147 \n "
    "PL P       15      30.974  1.8     1.06    255     165     0 \n "
    "NTL    N       7       14.007  1.55    0.75    0       191     255 \n "
    "CH     C       6       14.027  1.65    0.77    112     128     144 \n "
    "CH2    C       6       14.027  1.65    0.77    112     128     144 \n "
    "CH3    C       6       15.035  1.77    0.77    112     128     144 \n "
    "CE C       6       14.027  1.80    1.06    255     105     180  \n "
    "CK C       6       15.035  1.71    1.06    255     0     0  \n "
    "PO4      P       15      109.0  2.23     1.73    255     165     0 \n "
    "NC4     N       7       73.137  2.06    1.56    0       191     255 \n "
    "HDP       X       0       0.0     1.0     1.0     0     255     0  \n "
    "FAKE       X       0       0.0     1.0     1.0     0     255     0 \n "
    "SSD    X       0       0.0     1.0     1.0    255     20      147  \n "};

}

#endif

