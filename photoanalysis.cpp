#include "photoanalysis.h"
#include "ui_photoanalysis.h"
// #include <QImage>
#include <QDebug>

#include <algorithm>
#include "CImg.h"
//#include <iostream>

using namespace cimg_library;
//using namespace std;

template<typename T>
QList<QList<int>> maximas(const CImg<T> &t)
{
    QList<QList<int>>r;
    for (int z = 1; z < t.depth()-1; ++z){
        for (int y = 1; y < t.height()-1; ++y){
            for (int x = 1; x < t.width()-1; ++x){
                for (int c = 0; c < t.spectrum(); ++c){
                    QList<T> vals;
                    vals << t(x-1, y-1, z-1, c) << t(x-1, y, z-1, c) << t(x-1, y+1, z-1, c) << t(x, y-1, z-1, c) << t(x, y, z-1, c) << t(x, y+1, z-1, c) << t(x+1, y-1, z-1, c)
                         << t(x+1, y, z-1, c) << t(x+1, y+1, z-1, c)
                         << t(x-1, y-1, z, c) << t(x-1, y, z, c) << t(x-1, y+1, z, c) << t(x, y-1, z, c) << t(x, y+1, z, c) << t(x+1, y-1, z, c) << t(x+1, y, z, c)
                         << t(x+1, y+1, z, c)
                         << t(x-1, y-1, z+1, c) << t(x-1, y, z+1, c) << t(x-1, y+1, z+1, c) << t(x, y-1, z+1, c) << t(x, y, z+1, c) << t(x, y+1, z+1, c) << t(x+1, y-1, z+1, c)
                         << t(x+1, y, z+1, c) << t(x+1, y+1, z+1, c);
                    const auto val = t(x, y, z, c);
                    const auto min = *std::min_element(vals.begin(), vals.end());
                    const auto max = *std::max_element(vals.begin(), vals.end());
                    const T hysterisis = 0;//max-min;
                    if (((val+hysterisis*0.2) >= min && (val-hysterisis*0.2) <= max) || val != val)
                        continue;
                    QList<int>e;
                    e << x << y << z << val << c;
                    r << e;
                }
            }
        }
    }
    return r;
}

template<typename T>
QList<QList<T>> analyseGradient(const CImg<T> &t, const QString &nomFichier)
{
    auto base = t.get_RGBtoHSL();
    base.draw_image(0, 0, 0, 0, base.get_channel(0).normalize(0, 1));
    base.draw_image(0, 0, 0, 1, base.get_channel(1).normalize(0, 1));
    base.draw_image(0, 0, 0, 2, base.get_channel(2).normalize(0, 1));
    base.save((nomFichier + QString("base.bmp")).toLocal8Bit());
    auto grad = base.get_gradient();
    auto norm = (grad[0].get_sqr() + grad[1].get_sqr()).sqrt();
    auto phi = grad[1].get_atan2(grad[0]);
    phi.save((nomFichier + QString("phi.bmp")).toLocal8Bit());
    norm.save((nomFichier + QString("norm.bmp")).toLocal8Bit());

    QList<QList<T>>lPOIs;

    for (int z = 0; z < phi.depth(); ++z){
        for (int y = 0; y < phi.height(); ++y){
            for (int x = 0; x < phi.width(); ++x){
                // for (int c = 0; c < phi.spectrum(); ++c){
                if (norm(x, y, z, 0)<0.4) {
                    for (int c = 0; c < phi.spectrum(); ++c){
                        phi(x,y,z,c) = 0;
                    }
                }
                else
                {
                    QList<T>lPoi;
                    lPoi << x << y;

                    for (int c = 0; c < norm.spectrum(); ++c){
                        lPoi << norm(x, y, z, c);
                    }
                    for (int c = 0; c < phi.spectrum(); ++c){
                        lPoi << phi(x, y, z, c);
                    }
                    lPOIs << lPoi;
                }
                // }
            }
        }
    }
    phi.save((nomFichier + QString("filtered.bmp")).toLocal8Bit());

    return lPOIs;
}

PhotoAnalysis::PhotoAnalysis(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::PhotoAnalysis)
{
    ui->setupUi(this);
    //Les noms des fichiers qu'on va traiter
    QList<QString> fileNames;
    // Les noms des fichiers qu'on utilise
    //fileNames.push_back("C:\\Users\\sauvagee\\Documents\\GitHub\\SIFT\\imgs\\canape1.jpg");
    //fileNames.push_back("C:\\Users\\sauvagee\\Documents\\GitHub\\SIFT\\imgs\\canape2.jpg");
    // fileNames.push_back("/home/nsauvage/Images/Nausicaa.bmp");
    fileNames.push_back("kermit000.jpg");
    fileNames.push_back("kermit001.jpg");
    QList<QList<float>>prevPOIs;
    for (const auto & f : fileNames)
    {
        CImg<float> image(f.toLocal8Bit());
        auto name = f.split("/").last().split(".").first();
        auto lPOIs = analyseGradient(image, name);
        if (!prevPOIs.empty()) {
            testMatches(lPOIs, prevPOIs);
        }
        prevPOIs = lPOIs;
        // QList<QList<int>>pointsInteret;
        // for (auto ref = qMin(image.width(), image.height()), nbOctaves = 0;ref > 16; ref >>= 1, nbOctaves++) {
        //     CImg<float>slices;
        //     auto prev = image.get_RGBtoYCbCr().get_channel(0);
        //     for (int i = 0, k = 2; i < 10; ++i)
        //     {
        //         auto blurred = image.get_RGBtoYCbCr().get_channel(0);
        //         blurred = blurred.blur(2*pow(k, i));
        //         auto c = prev - blurred;
        //         blurred.save(QString("blur%1.bmp").arg(i).toLocal8Bit());
        //         prev.save(QString("prev%1.bmp").arg(i).toLocal8Bit());
        //         prev = blurred;
        //         slices = slices.append(c, 'z');
        //         auto out = c;
        //         out.save(QString("img%1Octave%2.bmp").arg(i).arg(nbOctaves).toLocal8Bit());
        //   }
        //     pointsInteret << maximas(slices);
        //     break;
        //     image = image.resize_halfXY();
        // }
        // image = CImg<float>(f.toLocal8Bit());
        // for (const auto &p: pointsInteret){
        //     image(p.at(0), p.at(1), 0, 0) = 255;
        //     image(p.at(0), p.at(1), 0, 1) = 0;
        //     image(p.at(0), p.at(1), 0, 2) = 255;
        // }
        // image.save("poi.bmp");
        ui->lblNbPoints->setText(QString("%1").arg(lPOIs.length()));
    }
}

PhotoAnalysis::~PhotoAnalysis()
{
    delete ui;
}

template<typename T>
void PhotoAnalysis::testMatches(const QList<QList<T> > &a, const QList<QList<T> > &b)
{
    for (const auto & pA : a)
    {
        auto gradA = pA;
        gradA.takeFirst();
        gradA.takeFirst();
        for (const auto & pB : b)
        {
            auto gradB = pB;
            gradB.takeFirst();
            gradB.takeFirst();
            bool equivalent = true;
            for (int i = 0; i < gradB.length(); ++i) {
                if (qAbs(gradB[i] - gradA[i]) > 0.1*qMax(qAbs(gradB[i]), qAbs(gradA[i]))) {
                    equivalent = false;
                    break;
                }
            }
            if (equivalent) {
                qDebug() << "Matching point : " << pA << pB;
            }
        }
    }
}
