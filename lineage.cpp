#include "lineage.h"
#include <QTextStream>
#include "simulation.h"
#include <QDebug>
#include <QList>
#include "genus.h"
#include "mainwindow.h"
#include <QSet>

/////////////////////////////////////////////////////////////////////
//Lineage class - equates essentially to a species. Lineages have
//a timestamp for when they appeared, went extinct, or split
//also holds pointers to parent and two daughter lineages, initialised to 0
//timestamps are initialised to -1, and remain so if unused (e.g. if lineage
//goes extinct and does not split, time_split remains -1. Split lineages
//are replaced by their daughters in the simulation
//dontpropogatedelete is a flag used in stripping extinct genera to stop daughter lineages being
//deleted before they can be reassigned.
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//Imports of globals from simulations.cpp
/////////////////////////////////////////////////////
extern double RDT_THRESHOLD;
extern int IDT_THRESHOLD;
extern int SCT_THRESHOLD;
extern int TCT_THRESHOLD;
extern int STT_THRESHOLD;
extern int EXTRA_MUTATES;
extern int FDTPLUS_THRESHOLD;
extern bool COUPLE_RATES;
extern double COUPLE_OFFSET;
extern double SPECIATION_MODIFIER;
extern double EXTINCTION_MODIFIER;
extern quint32 MUTATE_CHANCES[];
extern quint64 MUTATE_COUNTS[];
extern double fourier_b[];

extern qint64 genusnumberidt;
extern int character_mutation_mode;
extern int parameter_mode;
extern Lineage *dummy_parameter_lineage;

/////////////////////////////////////////////////////
//Constructor/destructor
/////////////////////////////////////////////////////

Lineage::Lineage(quint32 characters[], Lineage *parent, qint64 timestamp, quint32 initialchars[])
//initial characters is default 0 pointer, in which case none supplied, use current
{

    //tracked=false;
    //if (parent==0) tracked=true;
    time_created=timestamp;
    time_died=-1;
    time_split=-1;
    daughter_lineage_A=0;
    daughter_lineage_B=0;
    genusnumber=0;
    parent_lineage=parent;
    simple_id=-1;

    dontpropogatedelete=false;
    id=TheSimGlobal->nextid++;
    //copy characters
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        current_characters[i]=characters[i];
        if (initialchars)
            initial_characters[i]=initialchars[i];
        else
            initial_characters[i]=characters[i];
    }

    if (parent!=0)
    {
        speciate_prob=parent->speciate_prob;
        extinct_prob=parent->extinct_prob;
    }
    calculate_extinction_and_speciation();
}

Lineage::~Lineage()
{
    if (dontpropogatedelete) return; //switch to turn off deletion of children
    //Recursively delete all daughter lineages via their destructors;
    if (daughter_lineage_A!=0) delete daughter_lineage_A;
    if (daughter_lineage_B!=0) delete daughter_lineage_B;
}

/////////////////////////////////////////////////////
//Iteration - MBL algorithm implementation
/////////////////////////////////////////////////////

void Lineage::domutation()
{
    quint32 r,r2;

    bool mutated;
    //mutate characters here
    switch (character_mutation_mode)
    {
        case CHARACTER_MODE_SFM:
        if (TheSimGlobal->GetRandom16()<CHANCE_MUTATE)
        {
            r=TheSimGlobal->GetRandom(); //select a 32 bit word
            r&=(CHARACTER_WORDS-1);
            r2=TheSimGlobal->GetRandom();
            r2&=31; //select a bit
            current_characters[r] ^= tweakers[r2];
        }
        break;

        case CHARACTER_MODE_SGM:
        case CHARACTER_MODE_VGM:
        for (int i=0; i<CHARACTER_WORDS*32; i++)
            if (TheSimGlobal->GetRandom()<MUTATE_CHANCES[i])
            {
                MUTATE_COUNTS[i]++;
                current_characters[i/32] ^= tweakers[i%32];
            }
        break;

        case CHARACTER_MODE_EGM:
        mutated=false;
        for (int i=0; i<CHARACTER_WORDS*32; i++)
            if (TheSimGlobal->GetRandom()<PER_LINEAGE_MUTATE_CHANCES[i])
            {
                MUTATE_COUNTS[i]++;
                current_characters[i/32] ^= tweakers[i%32];
                mutated=true;
            }
        if (mutated) recalcmutatechances();

        break;

        default: break;
    }

    calculate_extinction_and_speciation();

}

void Lineage::random_walk_rates()
{
    //tweak extinction and speciation rates
    quint32 r=TheSimGlobal->GetRandom();
    r%=3;
    //r==0 means not move - 1/3 of the time
    if (r==1) //up
    {
        speciate_prob+=SPECIATION_CHANGE_PER_STEP;
        //cap at both envelope and 1
        if (speciate_prob>CHANCE_SPECIATE_DOUBLE+SPECIATION_MODIFIER) speciate_prob=CHANCE_SPECIATE_DOUBLE+SPECIATION_MODIFIER;
        if (speciate_prob>1.0) speciate_prob=1.0;
    }
    if (r==2) //down
    {
        speciate_prob-=SPECIATION_CHANGE_PER_STEP;
        //cap at both envelope and 1
        if (speciate_prob<CHANCE_SPECIATE_DOUBLE-SPECIATION_MODIFIER) speciate_prob=CHANCE_SPECIATE_DOUBLE-SPECIATION_MODIFIER;
        if (speciate_prob<0) speciate_prob=0;
    }


    if (COUPLE_RATES)
    {
        //just link extinction to speciation
        extinct_prob=qBound(0.0,speciate_prob-COUPLE_OFFSET,1.0);
    }
    else
    {
        //repeat for extinction
        r=TheSimGlobal->GetRandom();
        r%=3;
        //r==0 means not move - 1/3 of the time
        if (r==1) //up
        {
            extinct_prob+=EXTINCTION_CHANGE_PER_STEP;
            //cap at both envelope and 1
            if (extinct_prob>CHANCE_EXTINCT_DOUBLE+EXTINCTION_MODIFIER) extinct_prob=CHANCE_EXTINCT_DOUBLE+EXTINCTION_MODIFIER;
            if (extinct_prob>1.0) extinct_prob=1.0;
        }
        if (r==2) //down
        {
            extinct_prob-=EXTINCTION_CHANGE_PER_STEP;
            //cap at both envelope and 1
            if (extinct_prob<CHANCE_EXTINCT_DOUBLE-EXTINCTION_MODIFIER) extinct_prob=CHANCE_EXTINCT_DOUBLE-EXTINCTION_MODIFIER;
            if (extinct_prob<0) extinct_prob=0;
        }

    }


}

void Lineage::calculate_extinction_and_speciation()
//works out base chance of extinct and speciation from globals and genome
//non-genetic methods random walk the values
//For genetic methods - add or subtract according to bitcount of the genome. LS half of each word for speciate, MS for extinct
//Scale according to EXTINCTION_MODIFIER and SPECIATION_MODIFIER
//such that min is (base-modifier), max is (base+modifier). Clamp to 0-1 as well
{

    if (parameter_mode==PARAMETER_MODE_FIXED || parameter_mode==PARAMETER_MODE_GAMMA)
    {
        extinct_prob=CHANCE_EXTINCT_DOUBLE;
        speciate_prob=CHANCE_SPECIATE_DOUBLE;
        return;
    }

    if (parameter_mode==PARAMETER_MODE_LOCAL_NON_GENETIC)
    {
        random_walk_rates();
        return;
    }

    if (parameter_mode==PARAMETER_MODE_GLOBAL_NON_GENETIC)
    {
        if (dummy_parameter_lineage==0) return;
        speciate_prob=dummy_parameter_lineage->speciate_prob;
        extinct_prob=dummy_parameter_lineage->extinct_prob;
    }

    if (parameter_mode==PARAMETER_MODE_LOCAL)
    //genetic control of parameters
    {
        int bitcountspeciate=0;
        int bitcountextinct=0;
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            bitcountspeciate+=bitcounts[current_characters[i]%65536];
            bitcountextinct+=bitcounts[current_characters[i]/65536];
        }
        //totals will be in range 0-16*CHARACTER_WORDS
        speciate_prob=qBound(0.0,CHANCE_SPECIATE_DOUBLE+
            ((double)(8*CHARACTER_WORDS-bitcountspeciate))*SPECIATION_MODIFIER/(double)(CHARACTER_WORDS*8)
            ,1.0);

        if (COUPLE_RATES)
            extinct_prob=qBound(0.0,speciate_prob-COUPLE_OFFSET,1.0);
        else
        extinct_prob=qBound(0.0,CHANCE_EXTINCT_DOUBLE+
            ((double)(8*CHARACTER_WORDS-bitcountextinct))*EXTINCTION_MODIFIER/(double)(CHARACTER_WORDS*8)
            ,1.0);

        return;

    }

    if (parameter_mode==PARAMETER_MODE_GLOBAL)
    //as for local genetic control, but use the dummy lineage
    {
        if (dummy_parameter_lineage==0) return; //should only happen when constructing that lineage
        int bitcountspeciate=0;
        int bitcountextinct=0;
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            bitcountspeciate+=bitcounts[dummy_parameter_lineage->current_characters[i]%65536];
            bitcountextinct+=bitcounts[dummy_parameter_lineage->current_characters[i]/65536];
        }
        //totals will be in range 0-16*CHARACTER_WORDS
        speciate_prob=qBound(0.0,CHANCE_SPECIATE_DOUBLE+
            ((double)(8*CHARACTER_WORDS-bitcountspeciate))*SPECIATION_MODIFIER/(double)(CHARACTER_WORDS*8)
            ,1.0);

        if (COUPLE_RATES)
            extinct_prob=qBound(0.0,speciate_prob-COUPLE_OFFSET,1.0);
        else
        extinct_prob=qBound(0.0,CHANCE_EXTINCT_DOUBLE+
            ((double)(8*CHARACTER_WORDS-bitcountextinct))*EXTINCTION_MODIFIER/(double)(CHARACTER_WORDS*8)
            ,1.0);
        return;

    }
}

void Lineage::recalcmutatechances()
{
    //for EGM algorithm only
    double a[FOURIER_TERMS];
    double c[FOURIER_TERMS];

    //both are in range 0-32. Work out by excising parts of the genome - no need to bitcount, just take directly
    //We have an unknown number of 32 bit words. Need 5 bit chunks. Overlapping is fine
    //We insist that FOURIER_TERMS is no more than 6, so these all come from first word (6*5 = 30)
    //We also insist that CHARACTER_WORDS is at least 2, so two available words for this.
    quint32 andterm=32;
    quint32 word=current_characters[0];
    for (int i=0; i<FOURIER_TERMS; i++)
    {
        quint32 val=word&(andterm-1);
        val/=(andterm/32);
        a[i]=(double)val;
        andterm*=32;
    }

    andterm=32;
    word=current_characters[1];
    for (int i=0; i<FOURIER_TERMS; i++)
    {
        quint32 val=word&(andterm-1);
        val/=(andterm/32);
        c[i]=(double)val;
        andterm*=32;
    }

    TheSimGlobal->setFourierChances(&(a[0]),&(fourier_b[0]), &(c[0]),&(PER_LINEAGE_MUTATE_CHANCES[0]));

}

void Lineage::iterate(qint64 timestamp)
{
    //iterate a lineage - i.e. perform MBL calculations
    if (time_died!=-1) return; //it's dead - no need to do anything
    if (time_split!=-1) //it has split - iterate it's daughters instead
    {
        daughter_lineage_A->iterate(timestamp);
        daughter_lineage_B->iterate(timestamp);
        return;
    }

    domutation(); //also works out correct extinction/speciation rates

    //get a random double (0-1)
    quint32 r=TheSimGlobal->GetRandom();
    double rdouble=((double)r)/65536;
    rdouble/=65536;

    //if (id%100==0) qDebug()<<"CORE"<<timestamp<<id<<" ext "<<extinct_prob<< " spec "<<speciate_prob;
    //if (parent_lineage!=nullptr)
    //    if (parent_lineage->id%100==0) qDebug()<<"OFFSPRING OF CORE"<<timestamp<<id<<" ext "<<extinct_prob<< " spec "<<speciate_prob;

    if (rdouble<extinct_prob) //lineage went extinct
    {
        TheSimGlobal->leafcount--;
        time_died=timestamp; //record when this happened
        return;
    }

    if (rdouble<(speciate_prob+extinct_prob)) //lineage speciated
    {
        TheSimGlobal->leafcount++;
        time_split=timestamp; //record when this happened
        daughter_lineage_A=new Lineage(current_characters, this,timestamp); //create two new daughter lineages with same characters
        daughter_lineage_B=new Lineage(current_characters, this,timestamp);
        //if (tracked) daughter_lineage_A->tracked=true; //tracking into one split

        for (int i=0; i<EXTRA_MUTATES; i++) //extra mutations for PE model
        {
            daughter_lineage_A->domutation();
            daughter_lineage_B->domutation();
        }
        return;
    }
}


/////////////////////////////////////////////////////
//Recursive functions to report on/summarise tree
/////////////////////////////////////////////////////

int Lineage::count_alive()
{
    //recursively count how many descendents of this lineage are alive (don't have a time_died set) - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split - value is count_alive of both branches
        return daughter_lineage_A->count_alive()+daughter_lineage_B->count_alive();

    if (time_died!=-1) return 0;
    return 1;
}

int Lineage::count_extinct()
{
    //recursively count how many descendents of this lineage are extinct (DO have a time_died set) - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split
        return daughter_lineage_A->count_extinct()+daughter_lineage_B->count_extinct();

    if (time_died==-1) return 0; //should be impossible anyway - it can't be split and extinct
    return 1;
}

int Lineage::count_branched()
{
    //recursively count how many descendents of this lineage speciated - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split
        return 1+daughter_lineage_A->count_branched()+daughter_lineage_B->count_branched();

    return 0;
}

void Lineage::getextantlist(QList<Lineage *> *list)
{
    //appends all descending extant species to the list (including this one_
    if (time_split==-1 && time_died==-1) //extant
    {
        list->append(this);
    }
    else
    {
        if (time_split!=-1) //if it split - recurse onto daughters
        {
            daughter_lineage_A->getextantlist(list);
            daughter_lineage_B->getextantlist(list);
        }
    }
}

void Lineage::getleaflist(QLinkedList<Lineage *> *list) //includes extinct
{
    //appends all descending extant or extinct species to the list (including this one)
    if (time_split==-1) //a leaf
    {
        list->append(this);
    }
    else
    {
        daughter_lineage_A->getleaflist(list);
        daughter_lineage_B->getleaflist(list);
    }
}

bool Lineage::isThisADescendent(Lineage *desc)
{
    //recurseive function - does the lineage passed appear in descendents of this lineage?
    if (this==desc) return true;
    if (time_split==-1) return false; //no descendents and it wasn't this, so false
    else
    {
        bool b1=daughter_lineage_A->isThisADescendent(desc);
        bool b2=daughter_lineage_B->isThisADescendent(desc);
        if (b2 || b1) return true; else return false;
    }
}

bool Lineage::DoesThisLeafSetContainACladeDescendingFromMe(QSet<Lineage *> *potentialspecieslist)
{
    //er... you can guess this one
    if (time_split==-1)
    {
        if (potentialspecieslist->contains(this)) return true; //fine, I am in the clade
        else return false; //ah, no I'm not in the list, so I am an excluded species
    }
    else
    {
        bool b1=daughter_lineage_A->DoesThisLeafSetContainACladeDescendingFromMe(potentialspecieslist);
        bool b2=daughter_lineage_B->DoesThisLeafSetContainACladeDescendingFromMe(potentialspecieslist);
        if (b2 && b1) return true; else return false;
    }
}

bool Lineage::isThisAnAncestor(Lineage *desc)
{
    //recursive function - does the lineage passed appear in ancestors of this lineage?
    if (this==desc) return true;
    if (!parent_lineage) return false; //hit root - didn't find it
    return parent_lineage->isThisAnAncestor(desc);
}

void Lineage::resetGenusLabels()
{
    genusnumber=0;
    if (time_split!=-1)
    {
        daughter_lineage_A->resetGenusLabels();
        daughter_lineage_B->resetGenusLabels();
    }
}



/////////////////////////////////////////////////////
//Generate tree-file output recursively
/////////////////////////////////////////////////////

QString Lineage::tntstring()
{
    //recursively generate TNT-format text description of tree
    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("%ld",id);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        int A=daughter_lineage_A->count_alive();
        int B=daughter_lineage_B->count_alive();

        if (A==0 && B==0) s="";
        if (A>0 && B>0)
        {
            out<<"("<<daughter_lineage_A->tntstring()<<" "<<daughter_lineage_B->tntstring()<<")";
        }
        if (A>0 && B==0)
        {
            return daughter_lineage_A->tntstring();
        }
        if (B>0 && A==0)
        {
            return daughter_lineage_B->tntstring();
        }
        return s;
    }
}

QString Lineage::numbertree(int add)
{
    //recursively generate number-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        return QString("%1:%2").arg(simple_id+add).arg(bl);

    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"("<<daughter_lineage_A->numbertree(add)<<","<<daughter_lineage_B->numbertree(add)<<"):"<<bl;
        return s;
    }
}

QString Lineage::newickstring()
{
    //recursively generate Newick-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("g%ld_s%ld:%d",genusnumber,id,bl);
        else
            s.sprintf("{g%ld_s%ld}:%d",genusnumber,id,bl);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"("<<daughter_lineage_A->newickstring()<<","<<daughter_lineage_B->newickstring()<<"):"<<bl;
        return s;
    }
}

QString Lineage::xmlstring()
{
    //recursively generate xml-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("<clade><name>g%ld_s%ld</name><branch_length>%d</branch_length></clade>\n",genusnumber,id,bl);
        else
            s.sprintf("<clade><name>Ext_s%ld</name><branch_length>%d</branch_length></clade>\n",id,bl);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"<clade><name>n"<<id<<"</name><branch_length>"<<bl<<"</branch_length>"<<daughter_lineage_A->xmlstring()<<daughter_lineage_B->xmlstring()<<"</clade>\n";
        return s;
    }
}

void Lineage::dosimplenumbers()
{
    //for tree algorithms really. Provides secondary numbering system from 1-n, with no gaps
    if (time_split==-1) //a leaf
        simple_id=TheSimGlobal->nextsimpleIDnumber++;
    else
    {
        daughter_lineage_A->dosimplenumbers();
        daughter_lineage_B->dosimplenumbers();
    }
}

/////////////////////////////////////////////////////
//IDT and MIT algoritms
/////////////////////////////////////////////////////

int Lineage::nodelength(int currenttime)
{
    //branchlength to present or split or extinction
    if (time_split==-1)
    {
        if (time_died==-1)
            return currenttime-time_created;
        else
            return time_died-time_created;
    }
    else
        return time_split-time_created;
}

void Lineage::genuslabels_IDT(qint64 currentgenusnumber, int currenttime)
{
    // - look at the two lengths for daughters. if EITHER is over threshold, both get diffent genus numbers
    // - if this is not split, just add to node and return as before

    if (time_split==-1) // a leaf - extinct or otherwise - just add to current genus
    {
        genusnumber=currentgenusnumber;
        Genus *g;
        if (!(genera.contains(currentgenusnumber))) //genus doesn't exist - create and add to genera
        {
            g=new Genus;
            g->id=currentgenusnumber;
            genera.insert(g->id,g);
        }
        else
            g=genera.value(currentgenusnumber); //already exists

        g->species.append(this); //add this species to genus
    }
    else // a split.
    {
        int d1=daughter_lineage_A->nodelength(currenttime);
        int d2=daughter_lineage_B->nodelength(currenttime);

        if (d1>IDT_THRESHOLD || d2>IDT_THRESHOLD)
        {
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_A->genuslabels_IDT(currentgenusnumber,currenttime);
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_B->genuslabels_IDT(currentgenusnumber,currenttime);
        }
        else
        {
            //node was not too long - keep same genus number for both children
            daughter_lineage_A->genuslabels_IDT(currentgenusnumber,currenttime);
            daughter_lineage_B->genuslabels_IDT(currentgenusnumber,currenttime);
        }
    }
}



bool Lineage::nodelength_idtm(int currenttime)
{
    //is branchlength to present or split or extinction for this or daughers over the threshold?
    //if so, return false, otherwise return true;
    if (time_split==-1)
    {
        if (time_died==-1)
            {if ((currenttime-time_created)>IDT_THRESHOLD) return false; else return true;}
        else
            {if ((time_died-time_created)>IDT_THRESHOLD) return false; else return true;}
    }
    else //split
    {
        if ((time_split-time_created)>IDT_THRESHOLD) return false; //branch to offspring too long
        //check it's all OK from here on in?
        bool d1=daughter_lineage_A->nodelength_idtm(currenttime);
        if (d1==false) return false; //A wasn't OK, so no point looking at B
        return daughter_lineage_B->nodelength_idtm(currenttime); //A was OK, so return whatever B was
    }
}

void Lineage::genuslabels_IDTm(qint64 currentgenusnumber, int currenttime)
{
    // - look at the two lengths for daughters. if EITHER is over threshold OR any daughers are, both get diffent genus numbers
    // - if this is not split, just add to node and return as before

    if (time_split==-1) // a leaf - extinct or otherwise - just add to current genus
    {
        genusnumber=currentgenusnumber;
        Genus *g;
        if (!(genera.contains(currentgenusnumber))) //genus doesn't exist - create and add to genera
        {
            g=new Genus;
            g->id=currentgenusnumber;
            genera.insert(g->id,g);
        }
        else
            g=genera.value(currentgenusnumber); //already exists

        g->species.append(this); //add this species to genus
    }
    else // a split.
    {
        bool d1=daughter_lineage_A->nodelength_idtm(currenttime);
        bool d2=daughter_lineage_B->nodelength_idtm(currenttime);

        //
        if (d1==false || d2==false)
        {
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_A->genuslabels_IDTm(currentgenusnumber,currenttime);
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_B->genuslabels_IDTm(currentgenusnumber,currenttime);
        }
        else
        {
            //node was not too long - keep same genus number for both children
            daughter_lineage_A->genuslabels_IDT(currentgenusnumber,currenttime);
            daughter_lineage_B->genuslabels_IDT(currentgenusnumber,currenttime);
        }
    }
}



/////////////////////////////////////////////////////
//FDT and FDT+ algoritm
/////////////////////////////////////////////////////

void Lineage::genuslabels_distance(qint64 currentlabel, int treelength)
{
    //used to apply FDT to whole tree. Recursive, normally called on rootspecies first
    if (time_split==-1)
    {
        //Not split, so extant or extinct

        genusnumber=currentlabel; //put in a genus - applies also to extinct, though we do nothing with those

        if (time_died==-1) //extant
        {
            Genus *g;
            if (!(genera.contains(currentlabel))) //genus doesn't exist - create and add to genera
            {
                g=new Genus;
                g->id=currentlabel;
                genera.insert(g->id,g);
            }
            else
                g=genera.value(currentlabel); //already exists

            g->species.append(this); //add this species to genus
        }
        //do nothing if extinct
        return;
    }
    else
    {
        //it split - see if split is before or after the threshold
        if ((treelength-time_split) > ABS_THRESHOLD)
        {
            //before - so the daughters are in separate genera
            daughter_lineage_A->genuslabels_distance(TheSimGlobal->nextgenusnumber++,treelength);
            daughter_lineage_B->genuslabels_distance(TheSimGlobal->nextgenusnumber++,treelength);
        }
        else
        {
            //after - so in same genera. Set labels of all children
            daughter_lineage_A->setgenuslabels(currentlabel);
            daughter_lineage_B->setgenuslabels(currentlabel);
        }
    }
}

void Lineage::genuslabels_fdtplus(qint64 currentlabel, int treelength)
{
    //used to apply FDT to whole tree. Recursive, normally called on rootspecies first
    if (time_split==-1)
    {
        //Not split, so extant or extinct

        genusnumber=currentlabel; //put in a genus - applies also to extinct, though we do nothing with those

        if (time_died==-1) //extant
        {
            Genus *g;
            if (!(genera.contains(currentlabel))) //genus doesn't exist - create and add to genera
            {
                g=new Genus;
                g->id=currentlabel;
                genera.insert(g->id,g);
            }
            else
                g=genera.value(currentlabel); //already exists

            g->species.append(this); //add this species to genus
        }
        //do nothing if extinct
        return;
    }
    else
    {
        //it split - see if split is before or after the threshold - also split if either daughter length is over plus threshold
        int d=qMax(daughter_lineage_A->nodelength(treelength),daughter_lineage_B->nodelength(treelength));
        if ((treelength-time_split) > ABS_THRESHOLD || d>FDTPLUS_THRESHOLD)
        {
            //daughters are in separate genera
            daughter_lineage_A->genuslabels_fdtplus(TheSimGlobal->nextgenusnumber++,treelength);
            daughter_lineage_B->genuslabels_fdtplus(TheSimGlobal->nextgenusnumber++,treelength);
        }
        else
        {
            //after - so in same genera. Set labels of all children
            daughter_lineage_A->setgenuslabels(currentlabel);
            daughter_lineage_B->setgenuslabels(currentlabel);
        }
    }
}

void Lineage::setgenuslabels(qint64 currentlabel, bool labelfossils)
{
    //recursively set all species descending from here to the passed genus label
    genusnumber=currentlabel;
    if (time_split!=-1)
    {
        daughter_lineage_A->setgenuslabels(currentlabel,labelfossils);
        daughter_lineage_B->setgenuslabels(currentlabel,labelfossils);
    }
    else
    {
        if (time_died==-1 || labelfossils) //if labelfossils on - well label it even if it did die!
        {
            Genus *g;
            if (!(genera.contains(currentlabel)))
            {
                g=new Genus;
                g->id=currentlabel;
                genera.insert(g->id,g);
            }
            else
                g=genera.value(currentlabel);

            g->species.append(this);
        }
    }
}


/////////////////////////////////////////////////////
//Recursive functions to remove extinct taxa from tree prior to RDT. Returns pointer to new root.
/////////////////////////////////////////////////////

void Lineage::cull_dead_branches()
{
    //recurse through tree. Where we find a branch time but only a single branch - merge the two structures, deleting the upper one
    if (time_split!=-1) //do nothing if not a split
    {
        if (daughter_lineage_A==0 && daughter_lineage_B==0)
        {
            qDebug()<<"Internal error in cull_dead_branches. Oh dear.";
        }
        if (daughter_lineage_A==0) //must mean B was not 0 - both 0 case was culled in strip_extinct
        {
            time_died=daughter_lineage_B->time_died;
            time_split=daughter_lineage_B->time_split;
            Lineage *todelete=daughter_lineage_B;
            daughter_lineage_A=todelete->daughter_lineage_A;
            daughter_lineage_B=todelete->daughter_lineage_B;
            for (int i=0; i<CHARACTER_WORDS; i++) current_characters[i]=todelete->current_characters[i];
            if (daughter_lineage_A) daughter_lineage_A->parent_lineage=this;
            if (daughter_lineage_B) daughter_lineage_B->parent_lineage=this;

            todelete->dontpropogatedelete=true;
            delete todelete;
            //parent and time created stay the same
            cull_dead_branches(); //do this again on me
            return;
        }
        if (daughter_lineage_B==0) //must mean A was not 0 - both 0 case was culled in strip_extinct
        {
            time_died=daughter_lineage_A->time_died;
            time_split=daughter_lineage_A->time_split;
            Lineage *todelete=daughter_lineage_A;
            daughter_lineage_A=todelete->daughter_lineage_A;
            daughter_lineage_B=todelete->daughter_lineage_B;
            todelete->dontpropogatedelete=true;
            for (int i=0; i<CHARACTER_WORDS; i++) current_characters[i]=todelete->current_characters[i];
            if (daughter_lineage_A) daughter_lineage_A->parent_lineage=this;
            if (daughter_lineage_B) daughter_lineage_B->parent_lineage=this;
            delete todelete;
            //parent and time created stay the same
            cull_dead_branches(); //do this again on me
            return;
        }
        //if here - it's a normal split
        daughter_lineage_A->cull_dead_branches();
        daughter_lineage_B->cull_dead_branches();
    }
}

Lineage *Lineage::strip_extinct(Lineage *parent)
{
    //return 0 pointer if there are no extant descendants
    //otherwise return pointer to new Lineage structure with extant only descs
    if (time_died==-1 && time_split==-1) //extant, didn't split
    {
        //create new lineage for this, return;
        Lineage *l = new Lineage(current_characters,parent,time_created,initial_characters);
        l->id=id;
        //everything else can stay at default
        return l;
    }
    if (time_died!=-1) //extinct branch
    {
        return (Lineage *)0; //return null pointer - this branch is not real
    }
    if (time_split!=-1) //split node
    {
        Lineage *l = new Lineage(current_characters,parent,time_created,initial_characters);
        l->id=id;

        l->time_died=-1;
        l->time_split=time_split;
        l->daughter_lineage_A=daughter_lineage_A->strip_extinct(l);
        l->daughter_lineage_B=daughter_lineage_B->strip_extinct(l);
        //if both returned 0 - no extant descendants - we return 0 too
        if (l->daughter_lineage_A==0 && l->daughter_lineage_B==0)
        {
            delete l;
            return (Lineage *)0;
        }
//        if (l->daughter_lineage_A) l->daughter_lineage_A->parent_lineage=this;
//        if (l->daughter_lineage_B) l->daughter_lineage_B->parent_lineage=this;

        return l;
    }
    qDebug()<<"Error in strip extinct - unhandled case";
    return (Lineage *)0;

}


/////////////////////////////////////////////////////
//Functions used by RDT algoritm
/////////////////////////////////////////////////////


Lineage * Lineage::getsister()
{
    //determine the sister clade to this clade - returns 0 if no sister clade (root)
    if (parent_lineage==0) return 0;

    //must be one of parents daughters - find which
    if (parent_lineage->daughter_lineage_A==this)
        return parent_lineage->daughter_lineage_B;
    else
        return parent_lineage->daughter_lineage_A;
}

bool Lineage::RDT_check(quint64 currenttime)
{
    //checks to see whether this lineage could be incorporated into a genus
    //This function is documented in the manuscript
    if (time_split==-1) //already culled fossils, so must be a leaf
        return true;
    else
    {
        if ((int)(((double)(currenttime-time_split))/RDT_THRESHOLD)>=currenttime-time_created)
        {
            if (daughter_lineage_A->RDT_check(currenttime) && daughter_lineage_B->RDT_check(currenttime))
                return true;
            else
                return false;
        }
        else return false; //fails 50% age rule
    }
}

void Lineage::RDT_incorporate(Genus *g)
{
    //incorporate all tips descending from here into this genus
    if (time_split==-1)
    {
        g->species.append(this);
        genusnumber=g->id;
    }
    else
    {
        daughter_lineage_A->RDT_incorporate(g);
        daughter_lineage_B->RDT_incorporate(g);
    }
}

///////////////////////////////////////////////////////////////////////
//Functions used by Character-based taxonomy (SCT,TCT, STT) algorithms//
///////////////////////////////////////////////////////////////////////

void Lineage::genuslabels_SCTm(qint64 currentgenusnumber)
{
    // Recurse up tree
    // Each split - are both child nodes compatible with my genome? - recurse up and check.
    //If they are - recurse up and set genera. If they aren't - both child nodes become separate calls to main function with new genus codes. This will work with fossils included.

    if (time_split==-1) // a leaf - extinct or otherwise - just add to current genus
    {
        genusnumber=currentgenusnumber;
        Genus *g;
        if (!(genera.contains(currentgenusnumber))) //genus doesn't exist - create and add to genera
        {
            g=new Genus;
            g->id=currentgenusnumber;
            genera.insert(g->id,g);
        }
        else
            g=genera.value(currentgenusnumber); //already exists

        g->species.append(this); //add this species to genus
    }
    else // a split.
    {
        if (daughter_lineage_A->ismygenomecloseto(current_characters))
        if (daughter_lineage_B->ismygenomecloseto(current_characters))
        {
            //both all fine and dandy - so recurse up and set them all to this genus, and return
            daughter_lineage_A->setgenuslabels(currentgenusnumber,true);
            daughter_lineage_B->setgenuslabels(currentgenusnumber,true);
            return; //no need to go further up here
        }
        //if we get here - can't put them both in same genus. Split and continue up

        currentgenusnumber=++genusnumberidt;
        daughter_lineage_A->genuslabels_SCTm(currentgenusnumber);
        currentgenusnumber=++genusnumberidt;
        daughter_lineage_B->genuslabels_SCTm(currentgenusnumber);
    }
}

bool Lineage::ismygenomecloseto(quint32 comparecharacters[])
{

    if (time_split==-1) //a leaf
    {
        if (TheSimGlobal->distancebetween(comparecharacters,current_characters)>SCT_THRESHOLD) return false;
        else return true;
    }
    else
    {
        if (daughter_lineage_A->ismygenomecloseto(comparecharacters))
        if (daughter_lineage_B->ismygenomecloseto(comparecharacters))
            return true; //both OK
        return false; // at least one was not OK
    }
}

void Lineage::addtogenus(qint64 genusnum)
{
    genusnumber=genusnum;
    Genus *g;
    if (!(genera.contains(genusnum))) //genus doesn't exist - create and add to genera
    {
        g=new Genus;
        g->id=genusnum;
        genera.insert(g->id,g);
    }
    else
        g=genera.value(genusnum); //already exists

    g->species.append(this); //add this species to genus
}

/////////////////////////////////////////////////////
//Find clade of a certain size within the tree
/////////////////////////////////////////////////////

Lineage * Lineage::find_clade_with_precise_size(int preciseleafcount)
{
    if (time_split==-1) return 0; //no daughters
    int a=daughter_lineage_A->count_alive();
    if (a==preciseleafcount) return daughter_lineage_A;
    int b=daughter_lineage_A->count_alive();
    if (b==preciseleafcount) return daughter_lineage_B;
    if (a>preciseleafcount)
    {
        Lineage *r=daughter_lineage_A->find_clade_with_precise_size(preciseleafcount);
        if (r!=0) return r;
    }
    if (b>preciseleafcount)
    {
        Lineage *r=daughter_lineage_B->find_clade_with_precise_size(preciseleafcount);
        if (r!=0) return r;
    }

    //if it gets here it's not found big enough daughter clades
    return 0;
}

/////////////////////////////////////////////////////
//Debugging or reporting functions
/////////////////////////////////////////////////////

int Lineage::persistedfor()
{
    //returns number of steps the lineage persisted for until
    //simulation ended or it split or it died
    if (time_died!=-1)
        return time_died-time_created;
    if (time_split!=-1)
        return time_split-time_created;

    //still alive - can't report this, don't know total number of iterations here!
    return -1;

}

QString Lineage::dump()
{
    QString output;
    QTextStream out(&output);
    out<<"<br /><br />Lineage dump<br />";
    out<<"timstamp created: "<<time_created<<"<br />";
    if (parent_lineage==0) out<<"No parent lineage - root species<br />";
    else
        out<<"Parent lineage pointer: "<<parent_lineage<<"\n";

    if (time_died==-1 && time_split==-1) out<<"Status: Alive"<<"<br />";
    if (time_died!=-1) out<<"Status: Extinct, died at "<<time_died<<"<br />";
    if (time_split!=-1) out<<"Status: Split at "<<time_split<<")"<<"<br />";
    out<<"initial characters<br />";
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(initial_characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    out<<"<br />";
    out<<"current characters<br />";
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(current_characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    out<<"<br />";

    /*
    //DUMP CHILDREN MINMINAL
    if (time_split!=-1)
    {
        out<<"As it split, showing child characters too<br/>";
        out<<"Child A lasted for another "<<daughter_lineage_A->persistedfor()<<" iterations<br/>";
        out<<"initial characters CHILD A<br />";
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            QString binary;
            binary=QString::number(daughter_lineage_A->initial_characters[i],2);
            while (binary.length()<32)
                binary.prepend("0");
            out<<binary;
        }
        out<<"<br />";
        out<<"current characters  CHILD A<br />";
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            QString binary;
            binary=QString::number(daughter_lineage_A->current_characters[i],2);
            while (binary.length()<32)
                binary.prepend("0");
            out<<binary;
        }
        out<<"<br />";

        out<<"initial characters CHILD B<br />";
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            QString binary;
            binary=QString::number(daughter_lineage_B->initial_characters[i],2);
            while (binary.length()<32)
                binary.prepend("0");
            out<<binary;
        }
        out<<"<br />";

        out<<"Child B lasted for another "<<daughter_lineage_B->persistedfor()<<" iterations<br/>";
        out<<"current characters  CHILD B<br />";
        for (int i=0; i<CHARACTER_WORDS; i++)
        {
            QString binary;
            binary=QString::number(daughter_lineage_B->current_characters[i],2);
            while (binary.length()<32)
                binary.prepend("0");
            out<<binary;
        }
        out<<"<br />";

    }
    */
    return output;
}

QString Lineage::getcharactersasstring(quint32 characters[])
{
    QString output;
    QTextStream out(&output);
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    return output;
}

int Lineage::count_zeros(int word, quint32 mask)
{

    if (time_split!=-1) //split - recurse
    {
        return daughter_lineage_A->count_zeros(word,mask)+daughter_lineage_B->count_zeros(word,mask);
    }
    else
    {
        if ((current_characters[word]&mask)==0) return 1; else return 0;
    }
}

int Lineage::count_ones(int word, quint32 mask)
{

    if (time_split!=-1) //split - recurse
    {
        return daughter_lineage_A->count_ones(word,mask)+daughter_lineage_B->count_ones(word,mask);
    }
    else
    {
        if ((current_characters[word]&mask)==1) return 1; else return 0;
    }
}


QString Lineage::getcharactermatrix(bool root)
{
    QString retval;

    if (root)
    //this is root species - so outgroup - need starting characters
    {
        retval= QString("Species_%1        ").arg(0).left(14);
        retval+=getcharactersasstring(initial_characters);
        retval+="\n";
    }
    if (time_split!=-1) //split - recurse
    {
        return retval+daughter_lineage_A->getcharactermatrix(false)+daughter_lineage_B->getcharactermatrix(false);
    }
    return retval+QString("Species_%1        ").arg(simple_id).left(14)+getcharactersasstring(current_characters)+"\n";
}
