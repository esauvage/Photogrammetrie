#ifndef PHOTOANALYSIS_H
#define PHOTOANALYSIS_H

#include <QMainWindow>
//#include "CImg.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class PhotoAnalysis;
}
QT_END_NAMESPACE

class PhotoAnalysis : public QMainWindow
{
    Q_OBJECT

public:
    PhotoAnalysis(QWidget *parent = nullptr);
    ~PhotoAnalysis();

private:
    template<typename T>
    void testMatches(const QList<QList<T>> &a, const QList<QList<T>> &b);
    Ui::PhotoAnalysis *ui;
};
#endif // PHOTOANALYSIS_H
