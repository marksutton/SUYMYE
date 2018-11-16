#include "simulation.h"
#include "lineage.h"
#include "mainwindow.h"
#include <QString>
#include <QDebug>
#include <QTime>
#include <QFile>
#include <QMutableHashIterator>
#include <QMutableLinkedListIterator>
#include <QSet>
#include <QInputDialog>
#include "genus.h"
#include "math.h"
#include <random>

#define GAMMA_ALPHA 1
#define GAMMA_BETA .05

/////////////////////////////////////////////////////
//Simulation class - single-instance class to manage simulation
/////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//Global data declarations
/////////////////////////////////////////////////////

Simulation *TheSimGlobal;
QHash<qint64,Genus*> genera;
int CHANCE_EXTINCT;
int CHANCE_SPECIATE;
double CHANCE_SPECIATE_DOUBLE, CHANCE_EXTINCT_DOUBLE;

int ABS_THRESHOLD;
double RDT_THRESHOLD;
int IDT_THRESHOLD;
int CHANCE_MUTATE;
double dCHANCE_MUTATE;
double dCHANCE_MUTATE_MULTIPLIER;
int MCT_THRESHOLD;
int SCT_THRESHOLD;
int TCT_THRESHOLD;
int STT_THRESHOLD;
int FDTPLUS_THRESHOLD;

double SPECIATION_MODIFIER;
double EXTINCTION_MODIFIER;
double SPECIATION_CHANGE_PER_STEP;
double EXTINCTION_CHANGE_PER_STEP;


bool COUPLE_RATES;
double COUPLE_OFFSET;

int EXTRA_MUTATES;
quint32 MUTATE_CHANCES[CHARACTER_WORDS*32];
quint64 MUTATE_COUNTS[CHARACTER_WORDS*32];
double fourier_b[FOURIER_TERMS];
int character_mutation_mode;
int parameter_mode;

quint64 rtot;
quint64 rcount;
qint64 genusnumberidt;

quint32 tweakers[32]; //used for modifying characters during evolution
int bitcounts[65536]; //bitcount array - for each possible 16 bit word, how many bits are on?

Lineage *dummy_parameter_lineage;

/////////////////////////////////////////////////////
//Constructor/Destructor and run
/////////////////////////////////////////////////////

Simulation::Simulation()
{
    rootspecies=nullptr;
    crownroot=nullptr;
    TheSimGlobal=this;
    filepath="c:/";  //CHANGE THIS FOR NON-WINDOWS SYSTEMS!

    //set up list of data arrays - one for each possible type and a dummy for 0 (unclassified)
    for (int i=0; i<=TREE_MODE_MAX; i++) counts.append(new QHash<int,int>);
    for (int i=0; i<=TREE_MODE_MAX; i++) proportional_counts.append(new QList<int>);
    for (int i=0; i<=TREE_MODE_MAX; i++) saturated_tree_sizes.append(new QList<int>);
    for (int i=0; i<CHARACTER_WORDS*32; i++)  MUTATE_CHANCES[i]=0; //default - will get fixed at start of sim
    character_mutation_mode=CHARACTER_MODE_SFM;
    parameter_mode=PARAMETER_MODE_FIXED;
}

Simulation::~Simulation()
{
    //destructor - do nothing
}



void Simulation::run(MainWindow *mainwin)
{
    //performs a simulation run
    //keep pointer to main window instance
    mw=mainwin;


    std::default_random_engine generator;
    std::gamma_distribution<double> gd(GAMMA_ALPHA,GAMMA_BETA);

    stopflag=false;
    rfilepoint=0;

    //setup random stuff
    rtot=0;
    rcount=0;
    read65536numbers();

    rpoint=0;

    //do tweakers
    tweakers[0]=1;
    for (int n=1; n<32; n++) tweakers[n]=tweakers[n-1]*2;

    //and bitcounts array
    for (qint32 n=0; n<65536; n++)
    {
            qint32 count=0;
            for (int m=0; m<16; m++) if ((n & tweakers[m])!=0) ++count;  // count the bits
            bitcounts[n]=count;
    }

    //set up max genus size per tree
    QVector< QList<maxgenusdatapoint> > maxgenussize; //List of lists, done per taxonomy type. MS 32 bits - tree size. LS 32 bits - max genus size

    for (int i=0; i<=TREE_MODE_MAX; i++)
    {
        QList<maxgenusdatapoint> l;
        maxgenussize.append(l);
    }

    //clear any data from previous runs
    if (genera.count()>0)
        qDeleteAll(genera);
    genera.clear();
    for (int i=0; i<counts.length(); i++) counts.at(i)->clear(); //clear count data
    for (int i=0; i<proportional_counts.length(); i++)
    {
        proportional_counts.at(i)->clear();
        for (int j=0; j<=PROPORTIONAL_BINS; j++) //yes, <= is correct - one extra bin for 100%
            proportional_counts.at(i)->append(0); //set up empty bins
    }

    for (int i=0; i<saturated_tree_sizes.length(); i++) saturated_tree_sizes.at(i)->clear();
    //get settings from UI, in int and double types - int type is in range 0-65535
    CHANCE_EXTINCT=(int)(mw->getchanceextinct()*65535.0);
    CHANCE_SPECIATE=(int)(mw->getchancespeciate()*65535.0);

    CHANCE_EXTINCT_DOUBLE=mw->getchanceextinct();
    CHANCE_SPECIATE_DOUBLE=mw->getchancespeciate();
    SPECIATION_MODIFIER=mw->getspeciationmodifier();
    EXTINCTION_MODIFIER=mw->getextinctionmodifier();
    SPECIATION_CHANGE_PER_STEP=mw->getspeciationchangeperstep();
    EXTINCTION_CHANGE_PER_STEP=mw->getextinctionchangeperstep();

    COUPLE_RATES=mw->getCoupleRates();
    COUPLE_OFFSET=CHANCE_SPECIATE_DOUBLE-CHANCE_EXTINCT_DOUBLE;

    dCHANCE_MUTATE=mw->getchancemutate();
    dCHANCE_MUTATE_MULTIPLIER=mw->getmutatechancevariation();
    EXTRA_MUTATES=mw->getextramutations();
    STT_THRESHOLD=mw->getSTTthreshold();

    CHANCE_MUTATE=(int)(dCHANCE_MUTATE*65535.0);

    ABS_THRESHOLD=mw->getabsthreshold();
    RDT_THRESHOLD=mw->getRDTthreshold();
    IDT_THRESHOLD=mw->getIDTthreshold();
    SCT_THRESHOLD=mw->getSCTthreshold();
    TCT_THRESHOLD=mw->getTCTthreshold();
    FDTPLUS_THRESHOLD=mw->getFDTPLUSsplitthreshold();

    int preciseleafcount=mw->get_precise_leaf_count();

    double fourier_a[FOURIER_TERMS];
    double fourier_c[FOURIER_TERMS];
    for (int i=0; i<CHARACTER_WORDS*32; i++)  MUTATE_COUNTS[i]=0; //default - will get fixed at start of sim

    int generations=mw->getgenerations(); //this is number of trees to run
    int iterations=mw->getiterations(); //time iterations per tree
    int leafcountmax=mw->getmaxleafcount();
    bool usefossils=mw->getIncludeFossils();
    bool enforcemonophyly=mw->getEnforceMonophyly();
    parameter_mode=mw->getParameterMode();
    if ((parameter_mode==PARAMETER_MODE_LOCAL_NON_GENETIC || parameter_mode==PARAMETER_MODE_GLOBAL_NON_GENETIC || parameter_mode==PARAMETER_MODE_GLOBAL || parameter_mode==PARAMETER_MODE_LOCAL) && (SPECIATION_MODIFIER<0.00001 && EXTINCTION_MODIFIER<0.00001))
            parameter_mode=PARAMETER_MODE_FIXED; //if envelopes aren't set, just make it fixed

    QVector<quint32> saturation_array;
    saturation_array.fill(0,iterations*(CHARACTER_WORDS*32+1));


    //first up - refuse to run if certain menu combinations selected
    if (usefossils)
    {
        if (mw->getTaxonomyTypeInUse(TREE_MODE_FDT) || mw->getTaxonomyTypeInUse(TREE_MODE_FDT2) || mw->getTaxonomyTypeInUse(TREE_MODE_RDT))
        {
            mw->logtext("Error: taxonomy modes RDT, FDT and FDT+ are not compatible with the incorporation of extinct (fossil) leaves");
            return;
        }
    }
    mw->logtext("Starting simulation...");
    //run counts/totals
    int totalca=0; //count alive
    int totalce=0; //count extinct
    int totalcb=0; //count branched (speciated)
    int actualtreecount=0;
    int actualiterations=0;

    //main loop - iterate for specified number of trees
    for (int i=0; i<generations; i++)
    {
        actualiterations++;
        dummy_parameter_lineage=nullptr; //important - constructor will try to check it!

        if (parameter_mode==PARAMETER_MODE_GLOBAL)
        {
            quint32 dummycharacters[CHARACTER_WORDS];
            randomcharacters(dummycharacters); //set up characters for start
            //dummy lineage will be read for ext/speciation chance
            dummy_parameter_lineage=new Lineage(dummycharacters,(Lineage *)nullptr,0,nullptr); //run with single species that lacks a parent, created at time step 0
        }

        if (parameter_mode==PARAMETER_MODE_GLOBAL_NON_GENETIC)
        {
            quint32 dummycharacters[CHARACTER_WORDS];
            randomcharacters(dummycharacters); //set up characters for start
            //dummy lineage will be read for ext/speciation chance
            dummy_parameter_lineage=new Lineage(dummycharacters,(Lineage *)nullptr,0,nullptr); //run with single species that lacks a parent, created at time step 0

            dummy_parameter_lineage->extinct_prob=CHANCE_EXTINCT_DOUBLE;
            dummy_parameter_lineage->speciate_prob=CHANCE_SPECIATE_DOUBLE;
            //in this mode we are going to random walk extinction/speciation - so we set them to the specified start point

        }


        //sort out mutation mode - done for each tree.
        character_mutation_mode=mw->getCharacterMutationMode();
        if (character_mutation_mode==CHARACTER_MODE_SGM)
        {
            quint32 chance=(quint32)((double(65536)*(double)65536*mw->getchancemutate()));
            for (int i=0; i<CHARACTER_WORDS*32; i++)  MUTATE_CHANCES[i]=chance; //static over genome
        }

        if (character_mutation_mode==CHARACTER_MODE_VGM)
        {
            //set some random parameters. They are all in the range 0-31
            for (int i=0; i<FOURIER_TERMS; i++)
            {
                fourier_a[i]=(double)(GetRandom16()%32);
                fourier_c[i]=(double)(GetRandom16()%32);
                fourier_b[i]=0.5+4.0*((double)(GetRandom16())/65536.0);
            }
            setFourierChances(&(fourier_a[0]),&(fourier_b[0]),&(fourier_c[0]),&(MUTATE_CHANCES[0]));
        }

        if (character_mutation_mode==CHARACTER_MODE_EGM)
        {
            //set random parameters for fourier_b - these DON'T vary under genetic control, though a and c do
            for (int i=0; i<FOURIER_TERMS; i++)
            {
                fourier_b[i]=0.5+4.0*((double)(GetRandom16())/65536.0);
            }
        }

        if (mw->correct_number_trees_wanted()) i=0; //if iterating to correct number of trees ensure loop does not terminate on for statement

        //progress bar
        if (mw->correct_number_trees_wanted())
            mw->setProgress(actualtreecount,generations);
        else
            mw->setProgress(i,generations);

        //run data for run
        leafcount=1;
        nextgenusnumber=1;
        nextid=0;
        currenttime=0;
        if (rootspecies!=nullptr) delete rootspecies; //delete any existing root species from previous iteration of main loop
        if (crownroot!=nullptr) delete crownroot; //delete any existing fossil free tree from previous iteration of main loop
        crownroot=nullptr;
        //set up characters randomly
        quint32 startcharacters[CHARACTER_WORDS];
        randomcharacters(startcharacters);


        //create a root species
        rootspecies=new Lineage(startcharacters,(Lineage *)nullptr,0,nullptr); //run with single species that lacks a parent, created at time step 0
        rootspecies->extinct_prob=CHANCE_EXTINCT_DOUBLE;
        rootspecies->speciate_prob=CHANCE_SPECIATE_DOUBLE;


        //Run tree for correct number of iterations
        int j;


        for (j=0; j<iterations; j++)
        {
            if (parameter_mode==PARAMETER_MODE_GAMMA)
            {
                //select a new rate pair

                CHANCE_EXTINCT_DOUBLE=qBound(0.0,gd(generator),.495);
                if (COUPLE_RATES) CHANCE_SPECIATE_DOUBLE=CHANCE_EXTINCT_DOUBLE+0.01;
                else
                    CHANCE_SPECIATE_DOUBLE=qBound(0.0,gd(generator),.495)+0.01;
            }

            if (parameter_mode==PARAMETER_MODE_GLOBAL)
            {
                //do a mutation for the dummy parameter lineage
                dummy_parameter_lineage->domutation();
            }

            if (parameter_mode==PARAMETER_MODE_GLOBAL_NON_GENETIC)
            {
                //do a mutation for the dummy parameter lineage
                dummy_parameter_lineage->random_walk_rates();
                qDebug()<<j<<": speciation "<<dummy_parameter_lineage->speciate_prob<<"  extinction "<<dummy_parameter_lineage->extinct_prob;

            }

            increment();
            if (leafcount>leafcountmax) break;
        }

        if (dummy_parameter_lineage) delete dummy_parameter_lineage;

        //check - did it terminate because it hit leaf limit?
        if (j!=iterations)
        {
            //yes - log this
            QString so;
            QTextStream sout(&so);
            sout<<i<<": Leaf limit hit: Tree ran for "<<j<<" iterations, reached "<<leafcount<<" leaves.";
            mw->logtext(so,OUTPUT_VERBOSE);
            if (mw->throw_away_on_leaf_limit())
            {
                mw->logtext(QString("%1: Discarded (too many leaves)").arg(actualiterations), OUTPUT_VERBOSE);
                continue; //skip to next iteration of main loop
            }            
        }

        Lineage* useroot; //root to use for analyses - normally rootspecies, but not if preciseleafcount is set

        crownroot=rootspecies->strip_extinct((Lineage *)nullptr); //removes extinct tips
        if (crownroot!=nullptr)
        {

            crownroot->cull_dead_branches(); //merge single branch nodes

            if (preciseleafcount!=0) //cull to precise size if possible
            {

                Lineage *preciseleafcountroot = crownroot->find_clade_with_precise_size(preciseleafcount);
                if (preciseleafcountroot)
                {
                    //found one!
                    useroot=preciseleafcountroot;
                    crownroot=preciseleafcountroot; //both are the same precisely sizes clade
                }
                else
                {
                    mw->logtext(QString("%1: Discarded (no subtree of length %2)").arg(actualiterations).arg(preciseleafcount), OUTPUT_VERBOSE);
                    continue; //next iteration
                }
            }
            else
            {
                useroot=rootspecies;
            }

            if (mainwin->getSaturationNeeded()) DoSaturationData(useroot, j, &saturation_array);

            //get count of how many exant, extinct and branched lineages descend from the root species
            int ca=useroot->count_alive();
            int ce=useroot->count_extinct();
            int cb=useroot->count_branched();
            totalca+=ca; //add to totals
            totalce+=ce;
            totalcb+=cb;

            //log leaf count
            mw->logtext(QString("%1: extant leaves at end of run: %2").arg(actualiterations).arg(ca), OUTPUT_VERBOSE);


            //If there is a tree, i.e. we have any extant descendents of the root species

            actualtreecount++;

            //Produce version with extinct stripped

            nextsimpleIDnumber=0;
            crownroot->dosimplenumbers(); //do simple numbering scheme

            //Now go through the various taxonomy modes

            //1. TREE_MODE_UNCLASSIFIED
            mw->do_trees(TREE_MODE_UNCLASSIFIED,i,useroot); //handle the unclassified tree writing

            if (usefossils) leafcount=ca+ce; else leafcount=ca;

            //2. TREE_MODE_RDT
            if (mw->getTaxonomyTypeInUse(TREE_MODE_RDT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();
                //no need to clear classifications here - crownroot not yet used
                RDT_genera(); //do RDT classification

                int max=genera_data_report(TREE_MODE_RDT);
                maxgenusdatapoint dp;
                dp.maxgenussize=max;
                dp.treesize=leafcount;
                maxgenussize[TREE_MODE_RDT].append(dp);

                mw->do_trees(TREE_MODE_RDT,i,useroot);

                //No need to enforce monophyly - and not compatible with fossils - won't get this far with them turned on
            }

            //3. TREE_MODE_FDT
            if (mw->getTaxonomyTypeInUse(TREE_MODE_FDT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();
                crownroot->resetGenusLabels();
                crownroot->genuslabels_distance(0,j); //do FDT classification
                int max=genera_data_report(TREE_MODE_FDT);
                maxgenusdatapoint dp;
                dp.maxgenussize=max;
                dp.treesize=leafcount;
                maxgenussize[TREE_MODE_FDT].append(dp);
                mw->do_trees(TREE_MODE_FDT,i,useroot);
                //No need to enforce monophyly - and not compatible with fossils - won't get this far with them turned on
            }

            //3b. TREE_MODE_FDT+
            if (mw->getTaxonomyTypeInUse(TREE_MODE_FDT2))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();
                crownroot->resetGenusLabels();
                crownroot->genuslabels_fdtplus(0,j); //do FDT classification
                int max=genera_data_report(TREE_MODE_FDT2);
                maxgenusdatapoint dp;
                dp.maxgenussize=max;
                dp.treesize=leafcount;
                maxgenussize[TREE_MODE_FDT2].append(dp);
                mw->do_trees(TREE_MODE_FDT2,i,crownroot);
                //This doesn't use fossils currently, and no need to enforce monophyly
            }

            //4. TREE_MODE_IDT
            if (mw->getTaxonomyTypeInUse(TREE_MODE_IDT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();
                genusnumberidt=1;

                if (enforcemonophyly)
                {
                    //uses slightly different algorithm

                    if (usefossils) //use tree with extinct included
                    {
                        useroot->resetGenusLabels();
                        useroot->genuslabels_IDTm(genusnumberidt,j); //do IDT classification
                        int max=genera_data_report(TREE_MODE_IDT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_IDT].append(dp);
                        mw->do_trees(TREE_MODE_IDT,i,useroot);
                    }
                    else
                    {
                        crownroot->resetGenusLabels();
                        crownroot->genuslabels_IDTm(genusnumberidt,j); //do IDT classification
                        int max=genera_data_report(TREE_MODE_IDT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_IDT].append(dp);
                        mw->do_trees(TREE_MODE_IDT,i,crownroot);
                    }
                }
                else
                {
                    //without monophyly enforced
                    if (usefossils) //use tree with extinct included
                    {
                        useroot->resetGenusLabels();
                        useroot->genuslabels_IDT(genusnumberidt,j); //do IDT classification
                        int max=genera_data_report(TREE_MODE_IDT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_IDT].append(dp);
                        mw->do_trees(TREE_MODE_IDT,i,useroot);
                    }
                    else
                    {
                        crownroot->resetGenusLabels();
                        crownroot->genuslabels_IDT(genusnumberidt,j); //do IDT classification
                        int max=genera_data_report(TREE_MODE_IDT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_IDT].append(dp);
                        mw->do_trees(TREE_MODE_IDT,i,crownroot);
                    }
                }
            }

            //6. TREE_MODE_SCT //similarity characters
            if (mw->getTaxonomyTypeInUse(TREE_MODE_SCT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();

                if (enforcemonophyly)
                {
                    if (usefossils) //use tree with extinct included
                    {
                        useroot->resetGenusLabels();
                        useroot->genuslabels_SCTm(genusnumberidt); //do SCT mono classification
                        int max=genera_data_report(TREE_MODE_SCT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_SCT].append(dp);
                        mw->do_trees(TREE_MODE_SCT,i,useroot);
                    }
                    else
                    {
                        crownroot->resetGenusLabels();
                        crownroot->genuslabels_SCTm(genusnumberidt); //do SCT mono classification
                        int max=genera_data_report(TREE_MODE_SCT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_SCT].append(dp);
                        mw->do_trees(TREE_MODE_SCT,i,crownroot);
                    }
                }
                else //monophyly not enforced
                {
                    if (usefossils) //use tree with extinct included
                    {
                        useroot->resetGenusLabels();
                        do_SCT(useroot);
                        int max=genera_data_report(TREE_MODE_SCT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_SCT].append(dp);
                        mw->do_trees(TREE_MODE_SCT,i,useroot);
                    }
                    else
                    {
                        crownroot->resetGenusLabels();
                        do_SCT(crownroot);
                        int max=genera_data_report(TREE_MODE_SCT);
                        maxgenusdatapoint dp;
                        dp.maxgenussize=max;
                        dp.treesize=leafcount;
                        maxgenussize[TREE_MODE_SCT].append(dp);
                        mw->do_trees(TREE_MODE_SCT,i,crownroot);
                    }
                }
            }


            //7. TREE_MODE_TCT //type species algorithm
            if (mw->getTaxonomyTypeInUse(TREE_MODE_TCT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();

                if (usefossils) //use tree with extinct included
                {
                    useroot->resetGenusLabels();
                    do_TCT(useroot,enforcemonophyly);
                    int max=genera_data_report(TREE_MODE_TCT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_TCT].append(dp);
                    mw->do_trees(TREE_MODE_TCT,i,useroot);
                }
                else
                {
                    crownroot->resetGenusLabels();
                    do_TCT(crownroot,enforcemonophyly);
                    int max=genera_data_report(TREE_MODE_TCT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_TCT].append(dp);
                    mw->do_trees(TREE_MODE_TCT,i,crownroot);
                }
            }

            //9. TREE_MODE_STT //type species algorithm
            if (mw->getTaxonomyTypeInUse(TREE_MODE_STT))
            {
                //clear genus list
                if (genera.count()>0)  qDeleteAll(genera); //delete all existing genera data structures
                genera.clear();

                if (usefossils) //use tree with extinct included
                {
                    useroot->resetGenusLabels();
                    do_STT(useroot,enforcemonophyly);
                    int max=genera_data_report(TREE_MODE_STT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_STT].append(dp);
                    mw->do_trees(TREE_MODE_STT,i,useroot);
                }
                else
                {
                    crownroot->resetGenusLabels();
                    do_STT(crownroot,enforcemonophyly);
                    int max=genera_data_report(TREE_MODE_STT);
                    maxgenusdatapoint dp;
                    dp.maxgenussize=max;
                    dp.treesize=leafcount;
                    maxgenussize[TREE_MODE_STT].append(dp);
                    mw->do_trees(TREE_MODE_STT,i,crownroot);
                }
            }

            mw->plotcounts(&counts,false); //plot the results
        }
        else
            mw->logtext(QString("%1: Discarded (no extant leaves at end of simulation)").arg(actualiterations), OUTPUT_VERBOSE);
        //should we stop?
        if (mw->correct_number_trees_wanted() && actualtreecount==generations) stopflag=true;
        if (stopflag) {generations=i; break;}

    }

    mw->logtext("...Done!");

    //set progress bar to complete
    mw->setProgress(generations,generations);

    mw->plotcounts(&counts,true); //plot the results with basic frequency tables too - resets CSV

    mw->proportionaltables(&proportional_counts,&saturated_tree_sizes);  //this also starts the output file and writes freq tables, prop tables

    mw->outputmaxgenussizefile(&maxgenussize); //do the max genus stuff, appending to csv if needbe

    if (mainwin->getSaturationNeeded()) WriteSaturationData(iterations,&saturation_array,mainwin->getSaturationFileName());
}


/////////////////////////////////////////////////////
//Random number methods
/////////////////////////////////////////////////////

quint32 Simulation::GetRandom16()
{
    qint32 r=(quint32)((GetRandom() / 65536));
    return r;
}

quint32 Simulation::GetRandom()
{
    if (rpoint<65535)
    return randoms[rpoint++];
    else
    {
        read65536numbers();
        rpoint=0;
        return randoms[rpoint++];
    }
}

void Simulation::read65536numbers()
{
    //read more random numbers from file
    int c[10];
    rfilepoint+=65536*4;
#ifdef Q_OS_MAC  //should be outside package
    QFile f("../../../randomnumbers.dat");
#else //win32 or linux
    QFile f("randomnumbers.dat");
#endif

    //qDebug()<<"Random number file: "<<f.fileName();
    f.open(QIODevice::ReadOnly);
    if ((f.size()-rfilepoint)<65536*4) rfilepoint=qrand();
    f.seek(rfilepoint);
    f.read((char *)(&(randoms[0])),65536*4);
    f.close();

    for (int i=0; i<10; i++) c[i]=0;

    qint64 tot=0;
    for (int i=0; i<65536; i++)
    {
        double r=(double)randoms[i];
        r/=65536;
        r/=65536;
        r*=10;
        int ri=(int)r;
        c[ri]++;
        tot+=(qint64)randoms[i];
    }
}



/////////////////////////////////////////////////////
//Daughter methods for run
/////////////////////////////////////////////////////

void Simulation::increment()
{
    //one iteration of rootspecies
    currenttime++;
    if (rootspecies!=nullptr)
    {
        rootspecies->iterate(currenttime); //will chain-iterate them all
    }
}


void Simulation::do_SCT(Lineage *root)
{
    //Rule for Similarity Only
    //     - Gather all nodes (with or without fossils). Make list of all their genomes - this is the unassigned list
    //     - Compute distance matrix from each to each [no, do on fly]
    //     A. Pick a species at random, assign it to a genus, remove it from the unassigned list
    //          - Add all species with distance less than threshold into that genus, and remove them ALL from the list
    //          - For each of those species, add all ditto. Keep going until no more can be added
    //     - Repeat to A until list is empty.

    //1. Get list of leaves to classify
    qint64 thisgenusnumber=0;
    QLinkedList<Lineage *> leaves;
    root->getleaflist(&leaves);

    if (leaves.count()==0) return; //nothing - get out

    //set up data structures
    QList<Lineage *> *nextbatch=new QList<Lineage *>; //this holds the species moved to this genus this time
    QList<Lineage *> *thisbatch=new QList<Lineage *>;
    QList<Lineage *> *swaptemp;

    do //loop until all species are assigned
    {
        ++thisgenusnumber;
        //start with another seed - first one in the list
        Lineage *thisone=leaves.takeFirst();
        thisbatch->append(thisone);

        do
        {
            foreach (Lineage *l,(*thisbatch)) //for everything in this batch
            {
                l->addtogenus(thisgenusnumber);
                QMutableLinkedListIterator<Lineage *> it(leaves);
                while (it.hasNext())
                {
                    it.next();
                    Lineage *l2=it.value();
                    if (distancebetween(l->current_characters,l2->current_characters)<=SCT_THRESHOLD)
                    {
                        nextbatch->append(l2);
                        it.remove();
                    }
                };
            }
            swaptemp=thisbatch;
            thisbatch=nextbatch;
            nextbatch=swaptemp;
            nextbatch->clear();
        } while (thisbatch->count()!=0);
    } while (leaves.count()>0);

    //that OUGHT to be it!
    delete thisbatch;
    delete nextbatch;
}

void Simulation::do_STT(Lineage *root, bool enforcemonophyly)
{
    //Rules for Stratigraphic Type Species Taxonomy
    //exactly as TCT, but leaves also need to fall within a time-range of the type species
    //if fossils are excluded, this is the same as TCT

    qint64 thisgenusnumber=0;
    QLinkedList<Lineage *> leaves;
    root->getleaflist(&leaves);

    if (leaves.count()==0) return; //nothing - get out

    QMutableLinkedListIterator<Lineage *> it(leaves);
    //set up data structures
    do //loop until all species are assigned
    {
        ++thisgenusnumber;

        QSet<Lineage *> potentialspecieslist;

        //2. start with a randomly selected type species
        int pos=GetRandom()%leaves.count();
        it.toFront();
        it.next();
        for (int i=0; i<pos; i++) it.next();

        Lineage *typespecies=it.value();
        it.remove();

        potentialspecieslist.insert(typespecies);

        int tstime=typespecies->time_died;
        if (tstime==-1) tstime=currenttime;

        //3. Go through leaves, grab anything close enough and put in the genus
        it.toFront();
        while (it.hasNext())
        {
            it.next();
            Lineage *l2=it.value();

            if (distancebetween(typespecies->current_characters,l2->current_characters)<=TCT_THRESHOLD)
            {
                int cstime=l2->time_died;
                if (cstime==-1) cstime=currenttime;
                if (qAbs(cstime-tstime)<=STT_THRESHOLD)
                {
                    potentialspecieslist.insert(l2);
                    it.remove();
                }
            }
        };

        //4. [if enforce monophyly is on:] find largest clade within this genus that contains the type species
        //        - anything that doesn't fit in this clade - return to 'not placed' list
        if (enforcemonophyly)
        {
            //maintain a root of this clade node - start it at type species
            //iteratively:
            //  - look at parent of this node. If all children are in the genus, make this new root of this clade, continue
            //                                - if not, stop
            //Now remove all nodes that are NOT in this clade from the genus, and re-add to the 'to do' list
            Lineage *rootofclade=typespecies;
            bool done=false;
            do
            {
                done=true;
                if (rootofclade->parent_lineage)
                {
                    if (rootofclade->parent_lineage->DoesThisLeafSetContainACladeDescendingFromMe(&potentialspecieslist))
                    {
                        done=false;
                        rootofclade=rootofclade->parent_lineage;
                    }
                 }
            } while (!done);

            //remove species not in this clade
            foreach (Lineage *l, potentialspecieslist)
            {
                if (l->isThisAnAncestor(rootofclade))
                {
                   l->addtogenus(thisgenusnumber); //move it to new genus
                }
                else
                {
                    it.insert(l); //revert it to leaves list
                }
            }
        }
        else
        {
            //now actually assign the genus numbers
            foreach (Lineage *l, potentialspecieslist)
            {
                l->addtogenus(thisgenusnumber);
            }
        }
    } while (leaves.count()>0);
}


void Simulation::do_TCT(Lineage *root, bool enforcemonophyly)
{
    //Rules for Type Species Taxonomy

    //1. Gather all nodes into a list
    //2. Select an item at random from the list as type species, make a new genus
    //3. Grab things from the list that are within distance of this, place in genus
    //4. [if enforce monophyly is on:] find largest clade within this genus that contains the type species
    //        - anything that doesn't fit in this clade - return to 'not placed' list
    //5. Repeat 2-3 until nothing left in the list

    //1. Get list of leaves to classify
    qint64 thisgenusnumber=0;
    QLinkedList<Lineage *> leaves;
    root->getleaflist(&leaves);

    if (leaves.count()==0) return; //nothing - get out

    QMutableLinkedListIterator<Lineage *> it(leaves);
    //set up data structures
    do //loop until all species are assigned
    {
        ++thisgenusnumber;

        QSet<Lineage *> potentialspecieslist;

        //2. start with a randomly selected type species
        int pos=GetRandom()%leaves.count();
        it.toFront();
        it.next();
        for (int i=0; i<pos; i++) it.next();

        Lineage *typespecies=it.value();
        it.remove();

        potentialspecieslist.insert(typespecies);

        //3. Go through leaves, grab anything close enough and put in the genus
        it.toFront();
        while (it.hasNext())
        {
            it.next();
            Lineage *l2=it.value();
            if (distancebetween(typespecies->current_characters,l2->current_characters)<=TCT_THRESHOLD)
            {
                potentialspecieslist.insert(l2);
                it.remove();
            }
        };

        //4. [if enforce monophyly is on:] find largest clade within this genus that contains the type species
        //        - anything that doesn't fit in this clade - return to 'not placed' list
        if (enforcemonophyly)
        {
            //maintain a root of this clade node - start it at type species
            //iteratively:
            //  - look at parent of this node. If all children are in the genus, make this new root of this clade, continue
            //                                - if not, stop
            //Now remove all nodes that are NOT in this clade from the genus, and re-add to the 'to do' list
            Lineage *rootofclade=typespecies;
            bool done=false;
            do
            {
                done=true;
                if (rootofclade->parent_lineage)
                {
                    if (rootofclade->parent_lineage->DoesThisLeafSetContainACladeDescendingFromMe(&potentialspecieslist))
                    {
                        done=false;
                        rootofclade=rootofclade->parent_lineage;
                    }
                 }
            } while (!done);

            //remove species not in this clade
            foreach (Lineage *l, potentialspecieslist)
            {
                if (l->isThisAnAncestor(rootofclade))
                {
                   l->addtogenus(thisgenusnumber); //move it to new genus
                }
                else
                {
                    it.insert(l); //revert it to leaves list
                }
            }
        }
        else
        {
            //now actually assign the genus numbers
            foreach (Lineage *l, potentialspecieslist)
            {
                l->addtogenus(thisgenusnumber);
            }
        }
    } while (leaves.count()>0);
}

void Simulation::RDT_genera()
{
    if (crownroot==nullptr) return;
    //calculate genera using recursive RDT rules

    if (genera.count()>0) qDeleteAll(genera);
    genera.clear();

    QList<Lineage *>extantlist;
    crownroot->getextantlist(&extantlist);

    qint64 gnumber=1;

    for (int i=0; i<extantlist.count(); i++)
    {
        Lineage *thisone=extantlist[i];
        if (thisone->parent_lineage) //exclude root!
        {
            Lineage *sister=thisone->getsister();

            if (sister==nullptr) {continue;}

            if (sister->time_split==-1 && sister->genusnumber==0) //not split and must still be alive - so simple sister - though exclude if already labelled
            {
                //simple monophyletic pair. Create as a genus.
                Genus *g= new Genus;
                g->id=gnumber;
                g->species.append(thisone);
                g->species.append(sister);
                thisone->genusnumber=gnumber;
                sister->genusnumber=gnumber;
                g->rootnode=thisone->parent_lineage;
                gnumber++;
                genera.insert(g->id,g);

                bool continueloop=true;
                do
                {
                    //look at sister clade to this genus
                    sister=g->rootnode->getsister();


                    if (sister==nullptr)
                    {
                        //this happens if we try to expand out past root node. We can't - just move on
                        break;
                    }

                    int agegenusMRCA=currenttime-g->rootnode->time_split;
                    int agesistergenusMRCA=currenttime-g->rootnode->time_created;
                    if ((int)((double)agegenusMRCA/RDT_THRESHOLD)>=agesistergenusMRCA)
                    {
                        //passed the threshold cut off rule for depth - might be incorporatable
                        if (sister->RDT_check(currenttime))
                        {
                            sister->RDT_incorporate(g);
                            g->rootnode=g->rootnode->parent_lineage;
                            if (g->rootnode==crownroot) continueloop=false;
                        }
                        else
                            continueloop=false;
                    }
                    else
                        continueloop=false;
                }
                while (continueloop);
            }
        }
    }

    //deal with remaining singletons
    for (int i=0; i<extantlist.count(); i++)
    {
        if (extantlist[i]->genusnumber==0)
        {
            Genus *g= new Genus;
            g->id=gnumber;
            g->species.append(extantlist[i]);
            extantlist[i]->genusnumber=gnumber;
            gnumber++;
            genera.insert(g->id,g);
        }
    }


}

int Simulation::get_random_int_set(QSet<int> *extant_unused)
{
    int no=GetRandom()%extant_unused->count();

    int c=0;
    foreach(int i, *extant_unused)
    {
        if ((c++)==no) return i;
    }
    qDebug()<<"FELL OFF END OH DEAR";
    return -1;
}

int Simulation::get_mrca_age(QList<Lineage *> *terminals)
{
    //takes a list of nodes, and returns the age of the mrca.
    if (terminals->length()<2) qDebug("ERROR in get_mrca_age");
    Lineage *mrca=get_mrca(terminals->at(0),terminals->at(1));
    for (int i=2; i<terminals->length(); i++)
    {
        if (mrca->isThisADescendent(terminals->at(i))) continue; //skip
        else
            mrca=get_mrca(mrca,terminals->at(i));
    }
    return currenttime-mrca->time_split;
}

Lineage * Simulation::get_mrca(Lineage *l0, Lineage *l1)
{
    //determine most recent common ancestor (MRCA) for two lineages
    QList <Lineage *> l0ancestors, l1ancestors;

    while(l0->parent_lineage)
    {
        l0ancestors.append(l0->parent_lineage);
        l0=l0->parent_lineage;
    };
    while(l1->parent_lineage)
    {
        l1ancestors.append(l1->parent_lineage);
        l1=l1->parent_lineage;
    };

    //got lists, find first in L0 list that is also in L1 list

    foreach(Lineage *l, l0ancestors)
    {
        if (l1ancestors.contains(l)) return l;
    }

    //shouldn't get here!
    qDebug()<<"INTERNAL ERROR: No MRCA found in Simulation::get_mrca";
    return (Lineage *)nullptr;
}

int Simulation::genera_data_report(int mode)
//returns maximum genus size found
{

    int max=0;
    //add stuff to counts hash table for graphing and tables
    QHash<int,int> *thiscounts=counts[mode];
    foreach (Genus *g, genera)
    {
        int oldval=thiscounts->value(g->species.count(),0);
        (*thiscounts)[g->species.count()]=oldval+1;
    }

    QList<int> *thispcounts=proportional_counts[mode];
    foreach (Genus *g, genera)
    {
        if (g->species.count()>max) max=g->species.count();
        double proportion = ((double)g->species.count())/(double)leafcount;
        int bin=(int) (proportion * (100/PROPORTIONAL_BINS));
        if (g->species.count()==leafcount)
        {
            saturated_tree_sizes.at(mode)->append(leafcount); //store in list of trees where whole tree was one genus
            bin=PROPORTIONAL_BINS; //point to correct (final) bin for 100%
        }
        (*thispcounts)[bin]=thispcounts->at(bin)+1; //increment count in that bin
    }

    return max;
}

int Simulation::distancebetween(quint32 chars1[], quint32 chars2[])
//works out distance between a pair of genomes
{
    int total=0;
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        quint32 c1=chars1[i];
        quint32 c2=chars2[i];
        quint32 diffs = c1 ^ c2; //XOR the two to compare
        total+= bitcounts[diffs/(quint32)65536] +  bitcounts[diffs & (quint32)65535];
    }
    return total;
}

void Simulation::setFourierChances(double *a, double *b, double *c, quint32 *chance_array)
{
    //formula is Sum (a * sin ((b+n)*c*d)). a and c are 0-31. b is 0.5-5, precalced randoms for each run. d is PI()/(CHARATER_WORDS*2)
    double st=0;
    double dchances[CHARACTER_WORDS*32];
    double absmaxc=1.0;
    double d=M_PI/(CHARACTER_WORDS*2.0);
    //passed a set of fourier parameters, works out chance array - which is scaled to 32 bit uint

    //work out chancs as unscaled doubles
    for (int i=0; i<CHARACTER_WORDS*32; i++)
    {
        double ch=0;
        double di=(double)i;
        for (int j=0; j<FOURIER_TERMS; j++)
            ch+=a[j]*sin(b[j]*(di+c[j])*d);
        dchances[i]=ch;
        double ca=qAbs(ch);
        if (ca>absmaxc) absmaxc=ca; //track largest absolute value
    }

    //scale it all
    for (int i=0; i<CHARACTER_WORDS*32; i++)
    {
        double c=dchances[i];
        c/=absmaxc; //can't be 0
        c*=dCHANCE_MUTATE_MULTIPLIER;
        c+=dCHANCE_MUTATE;
        if (c<0) c=0;
        if (c>.999) c=.999;
        chance_array[i]=(quint32)(c*(65536.0)*(65536.0));
    }
}

/////////////////////////////////////////////////////
//Miscellaneous
/////////////////////////////////////////////////////

void Simulation::stop()
{
    stopflag=true; //triggers a stop at end of current tree
}


QString Simulation::modeToString(int mode)
{
    if (mode==TREE_MODE_FDT) return "Fixed Depth Taxonomy (FDT)";
    if (mode==TREE_MODE_RDT) return "Relative Top-Down Taxonomy (RDT)";
    if (mode==TREE_MODE_UNCLASSIFIED) return "No Taxonomy";
    if (mode==TREE_MODE_IDT) return "Internal Distance Taxonomy (IDT)";
    if (mode==TREE_MODE_FDT2) return "Fixed Depth Taxonomy Plus (FDT+)";
    if (mode==TREE_MODE_SCT) return "Similarity of Characters Taxonomy (SCT)";
    if (mode==TREE_MODE_TCT) return "Type-species Character Taxonomy (TCT)";
    if (mode==TREE_MODE_STT) return "Stratigraphic Type-species Taxonomy (STT)";
    return "error-typenotfound";
}

QString Simulation::modeToShortString(int mode)
{
    if (mode==TREE_MODE_FDT) return "FDT";
    if (mode==TREE_MODE_RDT) return "RDT";
    if (mode==TREE_MODE_UNCLASSIFIED) return "NoTax";
    if (mode==TREE_MODE_IDT) return "IDT";
    if (mode==TREE_MODE_FDT2) return "FDTPLUS";
    if (mode==TREE_MODE_SCT) return "SCT";
    if (mode==TREE_MODE_TCT) return "TCT";
    if (mode==TREE_MODE_STT) return "STT";
    return "error-typenotfound";
}

QString Simulation::dumpnewick(Lineage *rootlineage)
{
    return rootlineage->newickstring();
}

QString Simulation::dumptnt(Lineage *rootlineage)
{
    return "tread 'tree dumped from SUMYE in TNT format'\n"+rootlineage->tntstring()+"\nproc-;";
}

QString Simulation::dump_nex_alone(Lineage *rootlineage)
{
    QString ret;
    ret="#NEXUS\n";
    ret+="Begin trees;\n";
    ret+="    Translate\n";
    for (int i=0; i<nextsimpleIDnumber; i++)
    {
        QString comma;
        if (i!=(nextsimpleIDnumber-1)) comma=",";
        ret+=QString("%1 Species_%2%3\n").arg(i+1).arg(i).arg(comma);
    }
    ret+=";\n\n\n";

    QString tree=rootlineage->numbertree(1);
    //need to strip off last :number

    for (int i=tree.length()-1; i>0; i--)
    {
        if (tree[i]==QChar(':'))
        {
            tree=tree.left(i);
            break;
        }
    }
    ret+="tree tree1 =[&U]"+tree+";\n";
    ret+="END;\n";
    return ret;
}

QString Simulation::dumpphyloxml(Lineage *rootlineage)
{
    return "<?xml version=\"1.0\" encoding=\"UTF-8\"?><phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\"> <phylogeny rooted=\"true\" rerootable=\"true\">"
            + rootlineage->xmlstring()
            + "</phylogeny></phyloxml>\n";
}

void Simulation::randomcharacters(quint32 *chars)
//generates randomized characters for initial lineage
{
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        chars[i]=GetRandom();
    }
}



/////////////////////////////////////////////////////
//Functions for "saturation" comparison studies, looking at distance between taxa and genetic similarity
/////////////////////////////////////////////////////

void Simulation::WriteSaturationData(int iterations, QVector<quint32> *sat_array, QString fname)
{
    QFile satdata(fname);

    if (satdata.open(QIODevice::WriteOnly | QIODevice::Text)==false)
    {
        mw->logtext(QString("Couldn't open saturation file '%1' for output ").arg(fname));
        return;
    }
    QTextStream out(&satdata);

    //work out each depth line as 0-1 normalized value
    QVector<double> outputarray;
    outputarray.resize(iterations*(CHARACTER_WORDS*32+1));
    for (int i=0; i<iterations; i++)
    {
        quint64 tot=0;
        for (int j=0; j<CHARACTER_WORDS*32+1; j++)
            tot+=sat_array->at(i*(CHARACTER_WORDS*32+1)+j);

        //tot is now total occurrences of that particular depth of split (i)
        if (tot==0)  tot=1; //avoid divide by 0 errors - will be 0/tot anyway, so doesn't matter what value it has
        for (int j=0; j<CHARACTER_WORDS*32+1; j++)
            outputarray[i*(CHARACTER_WORDS*32+1)+j]=((double)sat_array->at(i*(CHARACTER_WORDS*32+1)+j))/((double)tot);

        //output array (same manual indexing as original) now has normalized 0-1 proportion of all possible similarities, i.e. a vertical line in output sums to 0
    }

    //output. Reverse order, from max to min similarity. x is depth, from 1 to (iterations), y is similarity, from 1 to 0 - there will be same numb
    for (int j=CHARACTER_WORDS*32; j>=0; j--)
    for (int i=0; i<iterations; i++)
    {
        if (i<iterations-1)
            out<<outputarray.at(i*(CHARACTER_WORDS*32+1)+j)<<",";
        else
            out<<outputarray.at(i*(CHARACTER_WORDS*32+1)+j)<<"\n";
    }

    out.flush();
    satdata.close();
}

void Simulation::DoSaturationData(Lineage *rootitem, int iterations, QVector<quint32> *sat_array)
{
    //First - assemble tip list

    QList <Lineage *> extantlist;
    rootitem->getextantlist(&extantlist);

    for (int i=0; i<(extantlist.count()-1); i++)
    {
        Lineage *l1=extantlist.at(i);
        for (int j=i+1; j<extantlist.count(); j++)  //pairwise comparison of all pairs in the extant list at end of sim
        {
            Lineage *l2=extantlist.at(j);

            Lineage *mrca=get_mrca(l1,l2);
            int depth=iterations-mrca->time_split;  //get age of split between the two


            //calculate genetic similarity between the two
            int similarity=CHARACTER_WORDS*32;
            for (int ii=0; ii<CHARACTER_WORDS; ii++)
            {
                quint32 diff=l1->current_characters[ii]^l2->current_characters[ii];
                similarity-=bitcounts[diff&65535];
                similarity-=bitcounts[diff/65536];
            }

            //Increment appropriate place in the array - which is actually a 2D array of depth of split and similarity, manually indexed (for speed I think)
            int index = depth*(CHARACTER_WORDS*32+1)+similarity;
            (*sat_array)[index]=sat_array->at(index)+1;

        }
    }
}
