#ifndef GENUS_H
#define GENUS_H
#include <QObject>
class Lineage;
class Genus
{
public:
    Genus();
    ~Genus();
    qint64 id;
    QList<Lineage *>species;
    Lineage *rootnode;
};

#endif // GENUS_H
