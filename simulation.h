#ifndef SIMULATION_H
#define SIMULATION_H
#include <QObject>
#include <QHash>
#include <QVector>

#define OUTPUT_QUIET 0
#define OUTPUT_VERBOSE 1

#define TREE_MODE_UNCLASSIFIED 0
#define TREE_MODE_FDT 1
#define TREE_MODE_RDT 2
#define TREE_MODE_FDT2 3
#define TREE_MODE_IDT 4
#define TREE_MODE_SCT 5
#define TREE_MODE_TCT 6
#define TREE_MODE_STT 7

#define CHARACTER_MODE_SFM 0
#define CHARACTER_MODE_SGM 1
#define CHARACTER_MODE_VGM 2
#define CHARACTER_MODE_EGM 3

#define PARAMETER_MODE_FIXED 0
#define PARAMETER_MODE_GAMMA 1
#define PARAMETER_MODE_GLOBAL 2
#define PARAMETER_MODE_LOCAL 3
#define PARAMETER_MODE_LOCAL_NON_GENETIC 4
#define PARAMETER_MODE_GLOBAL_NON_GENETIC 5

//Define below for TREE_MODE_MAX should be the max value used above
#define TREE_MODE_MAX 7

#define PROPORTIONAL_BINS 20 //100/this must be integral. So 10, 20, 25, 50 probably. Actually have this+1 bins - there is one for 100%

class Lineage;
class MainWindow;
class Genus;
class Simulation
{
public:
    Simulation();
    ~Simulation();
    void run(MainWindow *mainwin);
    MainWindow *mw;
    Lineage *rootspecies;
    Lineage *crownroot; //root for extant-only version of tree
    qint64 currenttime;
    void increment();
    quint32 randoms[65536];
    quint16 rpoint;
    quint64 rfilepoint;
    quint32 GetRandom();
    qint64 nextid;
    qint64 nextgenusnumber;
    bool stopflag;
    int leafcount;
    int nextsimpleIDnumber;

    void read65536numbers();
    quint32 GetRandom16();
    void stop();
    QString filepath;
    QString dumpnewick(Lineage *rootlineage);
    QString dumpphyloxml(Lineage *rootlineage);
    QString dumptnt(Lineage *rootlineage);
    QString modeToString(int mode);
    QList<QHash<int,int>*> counts;
    QList<QList<int>*> proportional_counts; //in PROPORTIONAL_BINS bins - 100/PROPORTIONAL_BINS must be an integer
    QList<QList<int>*> saturated_tree_sizes;
    QString modeToShortString(int mode);
    int distancebetween(quint32 chars1[], quint32 chars2[]);
    void setFourierChances(double *a, double *b, double *c, quint32 *chance_array);
    double gaussian_cdf(double x);
    double gaussian_pdf(double x);
    QString dump_nex_alone(Lineage *rootlineage);
private:
    quint32 rand256();

    int genera_data_report(int mode);

    void RDT_genera();
    Lineage *get_mrca(Lineage *l0, Lineage *l1);
    int get_random_int_set(QSet<int> *extant_unused);
    int get_mrca_age(QList<Lineage *> *terminals);
    void randomcharacters(quint32 *chars);
    void DoSaturationData(Lineage *rootitem, int iterations, QVector <quint32>*sat_array);
    void WriteSaturationData(int iterations, QVector <quint32>*sat_array, QString fname);
    void do_SCT(Lineage *root);
    void do_TCT(Lineage *root, bool enforcemonophyly);
    void do_STT(Lineage *root, bool enforcemonophyly);
};

struct maxgenusdatapoint
{
    quint32 maxgenussize;
    quint32 treesize;
};

extern Simulation *TheSimGlobal;
extern quint32 tweakers[32];
extern int bitcounts[65536];
extern QHash<qint64,Genus*> genera;
extern bool monomode;

#endif // SIMULATION_H
