// A little quick hack to do some rudimentry system profiling

#ifndef __MDPROFILE_H__
#define __MDPROFILE_H__

#define N_PROFILES 8
#define MAX_PROFILE_NAMELENGTH 40

enum proNames {pro1, pro2, pro3, pro4, pro5, pro6, pro7, pro8 };

extern void initProfile( void );

extern void startProfile( proNames theProfile );

extern void endProfile( proNames theProfile );

extern void writeProfiles( void );

#endif //__MDPROFILE_H__
