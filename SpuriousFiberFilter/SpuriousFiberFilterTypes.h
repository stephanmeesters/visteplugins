#ifndef bmia_SpuriousFiberFilterTypes_h
#define bmia_SpuriousFiberFilterTypes_h

namespace bmia {

/** Holding parameter settings */
typedef struct
{
    QString outputFiberDataName;

    double D33;
    double D44;
    double t;
    double windowSize;
    double epsilon;
    bool applyKernelInBothDirs;
    double minDist;
    bool requireRecompute;
    double cutoff;

    // saved results
    double avgScoreTotal;
    double* fiberMinScores;
    double* fiberScores;
    int* fiberStartIds;
    double goodness;

    bool applySubsampling;
    double samplingStep;

    bool applySelectAnterior;
    int numberOfAnteriorFibers;

} ParameterSettings;

}

#endif  // bmia_SpuriousFiberFilterTypes_h
