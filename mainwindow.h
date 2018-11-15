#ifndef MAINWINDOW_H
#define MAINWINDOW_H

class Lineage;
#include <QMainWindow>
#include <QHash>
#include <QActionGroup>
#include <QList>
#include <QVector>
#include <QStringList>
#include "simulation.h"

#define VERSION 3.0

namespace Ui {
class MainWindow;
}

class Simulation;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    Simulation *TheSim;

    bool gotsomedata;
    void logtext(QString text);
    void setProgress(int current, int max);
    double getchanceextinct();
    double getchancespeciate();
    double getchancemutate();
    double getthreshold();
    void plotcounts(QList<QHash<int, int> *> *data, bool showtable);
    int getgenerations();
    int getiterations();
    QActionGroup *alignmentGroup, *characterGroup, *parameterGroup;
    void do_trees(int mode, int iter, Lineage *rootlineage);
    int getmaxleafcount();
    int getabsthreshold();
    double getRDTthreshold();

    bool throw_away_on_leaf_limit();
    bool calculate_character_trees();
    bool per_tree_output();
    bool correct_number_trees_wanted();
    bool character_matrix();
    bool distance_matrix();
    int getIDTthreshold();
    int getSCTthreshold();
    int getMCTthreshold();
    int getTCTthreshold();
    bool getIncludeFossils();
    bool getEnforceMonophyly();

    bool getTaxonomyTypeInUse(int mode);
    int getCharacterMutationMode();
    double getmutatechancevariation();
    int getextramutations();
    int getSTTthreshold();
    int getFDTPLUSsplitthreshold();
    int getRTDPLUSmergethreshold();
    int getRTDPLUSsplitthreshold();
    double getextinctionmodifier();
    double getspeciationmodifier();
    void proportionaltables(QList<QList<int> *> *data, QList<QList<int> *> *satdata);
    bool getSaturationNeeded();
    QString getSaturationFileName();
    bool getTableFileNeeded();
    QString getCombinedTablesFileName();

    QList<QStringList> csvdata,csvdatabins;
    void outputmaxgenussizefile(QVector<QList<maxgenusdatapoint> > *data);
    void do_with_matrix_trees(int iter, Lineage *rootlineage);
    bool getCoupleRates();
    int getParameterMode();
    int get_precise_leaf_count();
    void log_identical_genomes(QString text);
public slots:
    void on_actionExit_triggered();
    void on_actionLogged_triggered();
    void on_actionStart_triggered();
    void on_actionStop_triggered();
    void on_actionSet_Export_Folder_2_triggered();
    void on_actionChart_to_PDF_triggered();
 private slots:
    void on_action_UseFDT_triggered();
    void on_action_UseIDT_triggered();
    void on_action_UseFDTPLUS_triggered();
    void on_action_UseRDTPLUS_triggered();
    void on_action_UseSCT_triggered();
    void on_action_UseRDT_triggered();
    void on_action_UseTCT_triggered();
    void on_action_UseSTT_triggered();
    void on_actionAbout_triggered();
    void on_action_UseUnclassified_triggered();

    void on_actionLog_scale_Y_axis_triggered();

private:
    Ui::MainWindow *ui;

    void setupGraphs();
    void addGraph(int index, QString colour, int treemode);
    QList<int> graphindices; //position is
    void outputTable(QString name, QVector<double> x, QVector<double> y, QHash<int, int> *count);
    void outputproportionalTable(QString name, QList<int> *data, QList<int> *satdata);
};

#endif // MAINWINDOW_H
