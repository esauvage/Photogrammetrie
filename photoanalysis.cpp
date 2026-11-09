#include "photoanalysis.h"
#include "ui_photoanalysis.h"
// #include <QImage>
#include <QDebug>
#include "CImg.h"

#include <algorithm>
#include "CImg.h"
//#include <iostream>

using namespace cimg_library;
//using namespace std;

float magenta[3] = {255, 0, 255};
float cyan[3] = {0, 255, 255};
float green[3] = {0, 255, 0};

template<typename T>
QList<float> PhotoAnalysis::descripteur(int x, int y, int z, int c, float angle, const CImg<T> &t) {
	QList <float> ret;
	CImg<float> vignette = t.get_slice(z).get_channel(c);
	// vignette.draw_arrow(x, y, x + 40 * cos(angle), y + 40 * sin(angle), magenta);
	// vignette.rotate(-angle*180./M_PI, x, y, 1);
	// vignette.save("vignette.bmp");
	vignette = t.get_slice(z).get_channel(c);
	vignette.rotate(-angle*180./M_PI, x, y, 1);
	auto dX = (vignette.width() - t.width()) / 2.;
	auto dY = (vignette.height() - t.height()) / 2.;
	auto L = vignette.get_gradient("x")[0];
	L = L.get_gradient("y")[0];
	for (int i = -6; i <= 6; i+=3) {
		for (int j = -6; j <= 6; j+=3) {
			ret << oriente(x + dX + i, y + dY + j, z, c, L, 8);
		}
	}
	float maxHisto = *std::max_element(ret.begin(), ret.end());
	for (auto & r : ret) {
		r /= maxHisto;
	}
	return ret;
}

template<typename T>
QList<float> PhotoAnalysis::oriente(int x, int y, int z, int c, const CImg<T> &L, int nbBins) {
	QList <float> ret;
	QList <float> bins;
	bins.assign(nbBins, 0);
	for (auto dX = -2; dX < 2; dX++) {
		for (auto dY = -2; dY < 2; dY++) {
		//Calcul de l'orientation
			auto Dxx1 = L.atXY(dX + x +1, dY + y, z, c);
			auto Dxx2 = L.atXY(dX + x -1, dY + y, z, c);
			auto Dyy1 = L.atXY(dX + x, dY + y + 1, z, c);
			auto Dyy2 = L.atXY(dX + x, dY + y - 1, z, c);
			double amplitude = sqrt((Dxx1-Dxx2)*(Dxx1-Dxx2) + ((Dyy1-Dyy2) * (Dyy1-Dyy2)));
			if (!amplitude) continue;
			double theta = atan((Dyy1-Dyy2)/(Dxx1-Dxx2));
			bins[theta/bins.size()] += amplitude;
		}
	}
	return bins;
}

template<typename T>
QList<QList<float>> PhotoAnalysis::maximas(const CImg<T> &DOFs)
{
	QList<QList<float>>r;
	for (int z = 1; z < DOFs.depth()-1; ++z){
		for (int y = 1; y < DOFs.height()-1; ++y){
			for (int x = 1; x < DOFs.width()-1; ++x){
				for (int c = 0; c < 1/*DOFs.spectrum()*/; ++c){
                    QList<T> vals;
					float dX = x;
					vals << DOFs(x-1, y-1, z-1, c) << DOFs(x-1, y, z-1, c) << DOFs(x-1, y+1, z-1, c) << DOFs(x, y-1, z-1, c) << DOFs(x, y, z-1, c) << DOFs(x, y+1, z-1, c) << DOFs(x+1, y-1, z-1, c)
						 << DOFs(x+1, y, z-1, c) << DOFs(x+1, y+1, z-1, c)
						 << DOFs(x-1, y-1, z, c) << DOFs(x-1, y, z, c) << DOFs(x-1, y+1, z, c) << DOFs(x, y-1, z, c) << DOFs(x, y+1, z, c) << DOFs(x+1, y-1, z, c) << DOFs(x+1, y, z, c)
						 << DOFs(x+1, y+1, z, c)
						 << DOFs(x-1, y-1, z+1, c) << DOFs(x-1, y, z+1, c) << DOFs(x-1, y+1, z+1, c) << DOFs(x, y-1, z+1, c) << DOFs(x, y, z+1, c) << DOFs(x, y+1, z+1, c) << DOFs(x+1, y-1, z+1, c)
						 << DOFs(x+1, y, z+1, c) << DOFs(x+1, y+1, z+1, c);
					const auto val = DOFs(x, y, z, c);
					const auto min = *std::min_element(vals.begin(), vals.end());
					const auto max = *std::max_element(vals.begin(), vals.end());
					const T hysterisis = 0;//max-min;
					if (((val+hysterisis*0.2) >= min && (val-hysterisis*0.2) <= max) || val != val)
						continue;
					auto D = DOFs.get_slice(z).get_gradient("xy");
					auto Dxx = D.at(0).get_gradient("x").at(0);
					auto Dyy = D.at(1).get_gradient("y").at(0);
					auto l1 = Dxx(x, y);
					auto l2 = Dyy(x, y);
					auto R = l1*l2 - 0.05*(l1 + l2)*(l1 + l2);//Criterion for corners
					if (R < 0.1)
						continue;
					// l1 = vals[10] - val;
					// l2 = vals[13] - val;
					// R = l1*l2 - 0.05*(l1 + l2)*(l1 + l2);//Criterion for corners
					// if (R < 0.1)
					// 	continue;
					// l1 = vals[15] - val;
					// l2 = vals[12] - val;
					// R = l1*l2 - 0.05*(l1 + l2)*(l1 + l2);//Criterion for corners
					// if (R < 0.1)
					// 	continue;
					// l1 = vals[10] - val;
					// l2 = vals[12] - val;
					// R = l1*l2 - 0.05*(l1 + l2)*(l1 + l2);//Criterion for corners
					// if (R < 0.1)
					// 	continue;
					double deltaX = 0.;
					do {
						auto d1 = DOFs.get_slice(z).get_channel(c).get_gradient("x").at(0);
						auto d2 = d1.get_gradient("x").at(0);
						d1.crop(dX, y-1, dX, y+1);
						d2.crop(dX, y-1, dX, y+1);
						d2.invert();
						// qDebug()<< d2(0, 0) << d2(0, 1) << d2(0, 2);
						// qDebug()<< d1(0, 0) << d1(1, 0) << d1(2, 0);
						deltaX = -(d2(0, 0) * d1(0, 0) + d2(0,1) * d1(1,0) + d2(0,2) * d1(2,0));
						dX += deltaX > 0.5 ? 1 : deltaX;
					} while (deltaX > 0.5);
					QList<float>e;
					e << dX << y << z << c << val;
					// auto save = DOFs.get_slice(0);
					// save.append(DOFs.get_slice(0), 'c');
					// save.append(DOFs.get_slice(0), 'c');
					// save.draw_circle(dX, y, 2, cyan);
					// save.save("Max.bmp");
					r << e;
                }
            }
        }
    }
    return r;
}

template<typename T>
QList<QList<T>> analyseGradient(const cimg_library::CImg<T> &t, const QString &nomFichier)
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
	// fileNames.push_back("C:\\Users\\sauvagee\\Documents\\GitHub\\SIFT\\imgs\\canape1.jpg");
	// fileNames.push_back("C:\\Users\\sauvagee\\Documents\\GitHub\\SIFT\\imgs\\canape2.jpg");
	// fileNames.push_back("/home/nsauvage/Images/Nausicaa.bmp");
	fileNames.push_back("kermit000.jpg");
	fileNames.push_back("kermit001.jpg");
	QList<QList<float>>prevPOIs;
	cimg_library::CImg<float> prevImage;
	for (const auto & f : fileNames)
    {
		cimg_library::CImg<float> image(f.toLocal8Bit());
		auto name = f.split("\\").last().split(".").first();
		// auto lPOIs = analyseGradient(image, name);
		QList<QList<float>>pointsInteret;
		for (auto ref = qMin(image.width(), image.height()), nbOctaves = 0;ref > 16; ref >>= 1, nbOctaves++) {
			cimg_library::CImg<float>DOGs;
			cimg_library::CImg<float>Ls;
			auto prev = image.get_RGBtoYCbCr().get_channel(0);
			for (int i = 0, k = 2; i < 10; ++i)
			{
				auto blurred = image.get_RGBtoYCbCr().get_channel(0);
				blurred = blurred.blur(2*pow(k, i));
				auto c = prev - blurred;
				// c.normalize(0, 255.);
				blurred.save(QString("blur%1.bmp").arg(i).toLocal8Bit());
				// prev.save(QString("prev%1.bmp").arg(i).toLocal8Bit());
				Ls = Ls.append(prev, 'z');
				Ls.get_slice(i).save(QString("prev%1.bmp").arg(i).toLocal8Bit());
				prev = blurred;
				DOGs = DOGs.append(c, 'z');
				DOGs.get_slice(i).normalize(0, 255.).save(QString("img%1Octave%2.bmp").arg(i).arg(nbOctaves).toLocal8Bit());
			}
			Ls = Ls.append(prev, 'z');
			auto POIs = maximas(DOGs);
			QList < QList <float>> toDel;
			for (auto& p : POIs) {
				p << nbOctaves;
				auto bins = oriente(p.at(0), p.at(1), p.at(2), p.at(3), Ls);
				auto angle = 0.;
				auto maxBin = 0;
				for (int i = 0; i < bins.size(); ++i) {
					if (bins[i] > maxBin) {
						angle = 360/bins.size() * (i + 0.5);
						maxBin = bins[i];
					}
				}
				p << angle;
				auto desc = descripteur(p.at(0), p.at(1), p.at(2), p.at(3), p.last(), Ls);
				p += desc;
				if (descMoy.isEmpty()) {
					descMoy = desc;
				} else {
					for (int i = 0; i < desc.size(); ++i) {
						descMoy[i] += desc[i];
					}
				}
			}
			for (auto &d : toDel) {
				POIs.removeAll(d);
			}
			pointsInteret << POIs;
			image = image.resize_halfXY();
		}
		for (int i = 0; i < descMoy.size(); ++i) {
			descMoy[i] /= pointsInteret.size();
		}
		double stdDev = 0;
		for (const auto &p: pointsInteret) {
			auto poi = p.last(descMoy.size());
			for (int i = 0; i < descMoy.size(); ++i) {
				stdDev += (poi[i] - descMoy[i])*(poi[i] - descMoy[i]);
			}
			stdDev = sqrt(stdDev);
		}
		// QList <QList <float>> toDel;
		// for (const auto &p1: pointsInteret) {
		// 	auto poi1 = p1.last(descMoy.size());
		// 	auto dev = 0.;
		// 	for (const auto &p2: pointsInteret) {
		// 		if (p1 == p2) continue;
		// 		auto poi2 = p2.last(descMoy.size());
		// 		for (int i = 0; i < descMoy.size(); ++i) {
		// 			dev += (poi1[i] - poi2[i])*(poi1[i] - poi2[i]);
		// 		}
		// 		dev = sqrt(dev);
		// 		if (dev < stdDev*0.1) {
		// 			toDel << p2 << p1;
		// 		}
		// 	}
		// }
		// qDebug() << "Descripteur moyen: " << descMoy;
		// for (auto& p : toDel) {
		// 	pointsInteret.removeAll(p);
		// }
		image = cimg_library::CImg<float>(f.toLocal8Bit());
		image.append(prevImage);
		for (const auto &p: pointsInteret){
			auto pX = p.at(0)*pow(2,p.at(5)) + 2 * p.at(6) * cos(p.at(7));
			pX = fmax(0, fmin(pX, image.width()));
			auto pY = p.at(1)*pow(2,p.at(5)) + 2 * p.at(6) * sin(p.at(7));
			pY = fmax(0, fmin(pY, image.height()));
			image.draw_arrow(p.at(0)*pow(2,p.at(5)), p.at(1)*pow(2,p.at(5)),
							 pX, pY, magenta);
			// for (int i = 9; i < p.size(); i+=1) {
			// 	pX = (p.at(0) + (i/8)/5)*pow(2,p.at(5));
			// 	pY = (p.at(1) + (i/8)%5)*pow(2,p.at(5));
			// 	image.draw_arrow(pX, pY, pX + 20 * cos(p.at(i)), pY + 20 * sin(p.at(i)), green);
			// }
		}
		ui->lblNbPoints->setText(ui->lblNbPoints->text() + QString("%1\n").arg(pointsInteret.length()));
		if (!prevPOIs.empty()) {
			auto matchs = testMatches(pointsInteret, prevPOIs);
			for (const auto &p: matchs) {
				image.draw_line(p.at(0)*pow(2,p.at(2)), p.at(1)*pow(2,p.at(2)),
								p.at(3)*pow(2,p.at(5))+ prevImage.width(), p.at(4)*pow(2,p.at(5)), cyan);
			}
		}
		image.save(QString("poi%1.bmp").arg(name).toLocal8Bit());
		ui->lblNbPoints->setPixmap(QPixmap(QString("poi%1.bmp").arg(name)));
		prevPOIs = pointsInteret;
		prevImage = image;
	}
}

PhotoAnalysis::~PhotoAnalysis()
{
	delete ui;
}

cimg_library::CImg<float> PhotoAnalysis::estimateF(const cimg_library::CImg<float> &A)
{
	return cimg_library::CImg<float>(3, 3, 1, 1).identity_matrix();
}

template<typename T>
QList<QList<T>> PhotoAnalysis::testMatches(const QList<QList<T> > &a, const QList<QList<T> > &b)
{
	auto wA = a;
	auto wB = b;
	cimg_library::CImg<float>fondamentalMatrix(9, 0, 1, 1);
    fondamentalMatrix.identity_matrix();
    QList<QList<T>> matches;
	for (int i = 0; i < fmin(fmin(16, a.size()), b.size()); ++i) {
		double minDist = INFINITY;
		QPair<int, int> bestPair;
		for (const auto & pA : wA)
		{
			const auto gradA = pA.last(52);
			for (const auto & pB : wB)
			{
				const auto gradB = pB.last(52);
				auto dist = 0.;
				for (int i = 0; i < gradB.length(); i+=2) {
					dist += (gradB[i+1] - gradA[i+1])*(gradB[i+1] - gradA[i+1]);
				}
				if (dist < minDist) {
					minDist = dist;
					bestPair.first = wA.indexOf(pA);
					bestPair.second = wB.indexOf(pB);
				}
			}
		}
		// if (equivalent) {
		qDebug() << "Matching point : " << minDist << "\n" << wA.at(bestPair.first)  << "\n" << wB.at(bestPair.second);
		QList<T> match;
		match << wA.at(bestPair.first)[0] << wA.at(bestPair.first)[1] << wA.at(bestPair.first)[5]
			  << wB.at(bestPair.second)[0] << wB.at(bestPair.second)[1] << wB.at(bestPair.second)[5];
		matches << match;
		wA.removeAll(wA.at(bestPair.first));
		wB.removeAll(wB.at(bestPair.second));
	}
	// cimg_library::CImg<float>row(9,1, 1, 1);
	// row(0,0) = match[2] * match[0];
	// row(1,0) = match[2] * match[1];
	// row(2,0) = match[2];
	// row(3,0) = match[3] * match[0];
	// row(4,0) = match[3] * match[1];
	// row(5,0) = match[3];
	// row(6,0) = match[0];
	// row(7,0) = match[1];
	// row(8,0) = 1;
	// fondamentalMatrix.append(row, 'y');
	// if (fondamentalMatrix.height() >= fondamentalMatrix.width()) {
	// 	auto solution = estimateF(fondamentalMatrix);
	// }
	return matches;
}
