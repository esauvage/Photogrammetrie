#ifndef PHOTOANALYSIS_H
#define PHOTOANALYSIS_H

#include <QMainWindow>
// #include "CImg.h"
namespace cimg_library {
template<typename T>
struct CImg;
}
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
	QList<QList<T>> testMatches(const QList<QList<T>> &a, const QList<QList<T>> &b);
	cimg_library::CImg<float> estimateF(const  cimg_library::CImg<float> &A);
	template<typename T>
	QList<QList<float>> maximas(const cimg_library::CImg<T> &DOFs);
	template<typename T>
	QList<float> oriente(int x, int y, int z, int c, const cimg_library::CImg<T> &t, int nbBins = 36);
	template<typename T>
	QList<float> descripteur(int x, int y, int z, int c, float angle, const cimg_library::CImg<T> &t);
	template<typename T>
	cimg_library::CImg<T> fondamentale(QList<QList<T> > matches);
	template<typename T>
	cimg_library::CImg<T> placePoint(const cimg_library::CImg<T> &P1, const cimg_library::CImg<T> &P2, cimg_library::CImg<T> x1, cimg_library::CImg<T> x2);
	void export2Ply(QList <QList <float>> points3D);

	Ui::PhotoAnalysis *ui;
	QList <float> descMoy;
};
#endif // PHOTOANALYSIS_H
