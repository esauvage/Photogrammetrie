#include "photoanalysis.h"
#include "ui_photoanalysis.h"
// #include <QImage>
#include <QDebug>
#include "CImg.h"

#include <algorithm>
#include <QFile>
#include <QRandomGenerator>
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
	for (int i = -8; i <= 8; i+=4) {
		for (int j = -8; j <= 8; j+=4) {
			QList <float> r;
			r.assign(8, 0);
			for (int l = -2; l <= 2; ++l) {
				for (int k = -2; k <= 2; ++k) {
					auto r11 = oriente(x + dX + i + l, y + dY + j + k, 0, 0, L, 8);
					for (int m = 0; m < r11.size(); ++m) {
						r[m] += r11[m];
					}
				}
			}
			ret << r;
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
	int degByBin = 360/nbBins;
	bins.assign(nbBins, 0);
	for (auto dX = -1; dX <= 1; dX++) {
		if ((dX + x + 1) >= L.width()) continue;
		if ((dX + x - 1) < 0) continue;
		for (auto dY = -1; dY <= 1; dY++) {
			if ((dY + y + 1) >= L.height()) continue;
			if ((dY + y - 1) < 0) continue;
			//Calcul de l'orientation
			auto Dxx1 = L(dX + x +1, dY + y, z, c);
			auto Dxx2 = L(dX + x -1, dY + y, z, c);
			auto Dyy1 = L(dX + x, dY + y + 1, z, c);
			auto Dyy2 = L(dX + x, dY + y - 1, z, c);
			double amplitude = sqrt((Dxx1-Dxx2)*(Dxx1-Dxx2) + ((Dyy1-Dyy2) * (Dyy1-Dyy2)));
			if (!amplitude) continue;
			int theta = atan((Dyy1-Dyy2)/(Dxx1-Dxx2)) * 180/M_PI + 360;
			theta %= 360;
			bins[theta/degByBin] += amplitude;
		}
	}
	// auto vignette = L.get_slice(z);
	// vignette.draw_circle(x, y, 2, magenta);
	// vignette.save("vignette.bmp");
	return bins;
}

template<typename T>
QList<QList<float>> PhotoAnalysis::maximas(const CImg<T> &DOFs)
{
	QList<QList<float>>r;
	for (int z = 1; z < DOFs.depth()-1; ++z){
		auto D = DOFs.get_slice(z).get_gradient("xy");
		auto Dxx = D.at(0).get_gradient("x").at(0);
		auto Dyy = D.at(1).get_gradient("y").at(0);
		for (int y = 1; y < DOFs.height()-1; ++y){
			for (int x = 1; x < DOFs.width()-1; ++x){
				for (int c = 0; c < 1/*DOFs.spectrum()*/; ++c){
					// QList<T> vals;
					auto vals = DOFs.get_crop(x-1,y-1, z-1, x+1, y+1, z+1);
					float dX = x;
					// vals << DOFs(x-1, y-1, z-1, c) << DOFs(x-1, y, z-1, c) << DOFs(x-1, y+1, z-1, c) << DOFs(x, y-1, z-1, c) << DOFs(x, y, z-1, c) << DOFs(x, y+1, z-1, c) << DOFs(x+1, y-1, z-1, c)
					// 	 << DOFs(x+1, y, z-1, c) << DOFs(x+1, y+1, z-1, c)
					// 	 << DOFs(x-1, y-1, z, c) << DOFs(x-1, y, z, c) << DOFs(x-1, y+1, z, c) << DOFs(x, y-1, z, c) << DOFs(x, y+1, z, c) << DOFs(x+1, y-1, z, c) << DOFs(x+1, y, z, c)
					// 	 << DOFs(x+1, y+1, z, c)
					// 	 << DOFs(x-1, y-1, z+1, c) << DOFs(x-1, y, z+1, c) << DOFs(x-1, y+1, z+1, c) << DOFs(x, y-1, z+1, c) << DOFs(x, y, z+1, c) << DOFs(x, y+1, z+1, c) << DOFs(x+1, y-1, z+1, c)
					// 	 << DOFs(x+1, y, z+1, c) << DOFs(x+1, y+1, z+1, c);
					const auto val = DOFs(x, y, z, c);
					// const auto min = *std::min_element(vals.begin(), vals.end());
					// const auto max = *std::max_element(vals.begin(), vals.end());
					const auto min = vals.min();
					const auto max = vals.max();
					const T hysterisis = 0;//max-min;
					if (((val+hysterisis*0.2) > min && (val-hysterisis*0.2) < max) || val != val)
						continue;
					auto l1 = Dxx(x, y);
					auto l2 = Dyy(x, y);
					auto R = l1*l2 - 0.05*(l1 + l2)*(l1 + l2);//Criterion for corners
					if (R < 0.1)
						continue;
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
			CImg<float>DOGs;
			CImg<float>Ls;
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
			for (auto& p : POIs) {
				p << nbOctaves;
				auto bins = oriente(p.at(0), p.at(1), p.at(2), p.at(3), Ls);
				auto angle = 0.;
				auto maxBin = 0;
				for (int i = 0; i < bins.size(); ++i) {
					if (bins[i] > maxBin) {
						angle = (2 * M_PI)/bins.size() * (i + 0.5);
						maxBin = bins[i];
					}
				}
				p << angle;
				auto desc = descripteur(p.at(0), p.at(1), p.at(2), p.at(3), p.last(), Ls);
				p += desc;
				// if (descMoy.isEmpty()) {
				// 	descMoy = desc;
				// } else {
				// 	for (int i = 0; i < desc.size(); ++i) {
				// 		descMoy[i] += desc[i];
				// 	}
				// }
			}
			pointsInteret << POIs;
			image = image.resize_halfXY();
		}
		// for (int i = 0; i < descMoy.size(); ++i) {
		// 	descMoy[i] /= pointsInteret.size();
		// }
		// double stdDev = 0;
		// for (const auto &p: pointsInteret) {
		// 	auto poi = p.last(descMoy.size());
		// 	for (int i = 0; i < descMoy.size(); ++i) {
		// 		stdDev += (poi[i] - descMoy[i])*(poi[i] - descMoy[i]);
		// 	}
		// 	stdDev = sqrt(stdDev);
		// }
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
			auto pX = p.at(0)*pow(2,p.at(5)) + 20 * cos(p.at(6));
			pX = fmax(0, fmin(pX, image.width()));
			auto pY = p.at(1)*pow(2,p.at(5)) + 20 * sin(p.at(6));
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
			CImg<float> definitiveF;
			float distMin = INFINITY;
			CImg<float>X1(2, 1, 1, 1);
			CImg<float>X2(2, 1, 1, 1);
			CImg<float> x1, x2, X;
			CImg<float> F;
			CImg<float> U, S, V;
			CImg<float> e;
			CImg<float> ep, P1, Proj2;
			CImg<float> skew(3,3,1,1,0);
			for (int i = 0; i < 100000; ++i) {
				QList<QList <float>> inliers_sample;
				std::sample(matchs.cbegin(), matchs.cend(),
							std::back_inserter(inliers_sample), 8,
							*QRandomGenerator::global());
				F = fondamentale(inliers_sample);
				F.SVD(U, S, V);
				// e  : vecteur propre droit de F  -> dernière colonne de V
				e  = V.get_column(2);
				// e' : vecteur propre gauche de F -> dernière colonne de U
				ep = U.get_column(2);
				// Normalisation (le signe est arbitraire)
				e  /= e(2);
				ep /= ep(2);
				// qDebug() << "||F*e||   = " << (F * e).magnitude(); // ~1e-15
				// qDebug() << "||F^T*e'||= " << (F.get_transpose() * ep).magnitude(); // ~1e-15
				P1 = F;
				P1.identity_matrix();
				P1.crop(0, 0, 3, 2);
				// qDebug() << "P1 :";
				// cimg_forX(P1, x) {
				// 	qDebug() << P1(x, 0) << P1(x, 1) << P1(x, 2) << P1(x, 3);
				// }

				skew(0,1) = -ep(2); skew(0,2) =  ep(1);
				skew(1,0) =  ep(2); skew(1,2) = -ep(0);
				skew(2,0) = -ep(1); skew(2,1) =  ep(0);

				// P2 = F;
				Proj2 = skew * F;
				Proj2.append(ep, 'x');
				Proj2(0, 3) = ep(0);
				Proj2(1, 3) = ep(1);
				Proj2(2, 3) = ep(2);
				// qDebug() << "P2 :";
				// cimg_forX(P2, x) {
				// 	qDebug() << P2(x, 0) << P2(x, 1) << P2(x, 2) << P2(x, 3);
				// }

				QList <QList <float>> points3D;
				float dist = 0;
				for (const auto &p: inliers_sample) {
					// QList<float>X1;
					X1(0) = p.at(0)*pow(2,p.at(2));
					X1(1) = p.at(1)*pow(2,p.at(2));
					X2(0) = p.at(3)*pow(2,p.at(5));
					X2(1) = p.at(4)*pow(2,p.at(5));
					// X1 << p.at(0)*pow(2,p.at(2)) << p.at(1)*pow(2,p.at(2));
					// QList<float>X2;
					// X2 << p.at(3)*pow(2,p.at(5)) << p.at(4)*pow(2,p.at(5));
					X = placePoint(P1, Proj2, X1, X2);
					qDebug() << X(0) << X(1) << X(2);
					x1 = P1*X;
					x2 = Proj2*X;
					dist += sqrt((X1(0)-x1(0)) * (X1(0)-x1(0)) + (X1(1)-x1(1)) * (X1(1)-x1(1)) + (X2(0)-x2(0)) * (X2(0)-x2(0)) + (X2(1)-x2(1)) * (X2(1)-x2(1)));
					QList <float> p3D;
					p3D << X(0) << X(1) << X(2);
					points3D << p3D;
				}
				if (dist < distMin) {
					distMin = dist;
					definitiveF = F;
					export2Ply(points3D);
				}
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
    QList<QList<T>> matches;
	for (int i = 0; wA.size() && wB.size() && matches.size()<100; ++i) {
		double minDist = INFINITY;
		double minDist2 = INFINITY;
		QPair<int, int> bestPair;
		QPair<int, int> bestPair2;
		for (const auto & pA : wA)
		{
			const auto gradA = pA.last(8*25);
			for (const auto & pB : wB)
			{
				const auto gradB = pB.last(8*25);
				auto dist = 0.;
				for (int i = 0; i < gradB.length(); i+=1) {
					dist += (gradB[i] - gradA[i])*(gradB[i] - gradA[i]);
				}
				if (dist < minDist) {
					minDist2 = minDist;
					minDist = dist;
					bestPair2 = bestPair;
					bestPair.first = wA.indexOf(pA);
					bestPair.second = wB.indexOf(pB);
				}
			}
		}
		QList<T> match;
		match << wA.at(bestPair.first)[0] << wA.at(bestPair.first)[1] << wA.at(bestPair.first)[5]
			  << wB.at(bestPair.second)[0] << wB.at(bestPair.second)[1] << wB.at(bestPair.second)[5];
		wA.removeAll(wA.at(bestPair.first));
		if (minDist/minDist2 > 0.99) {
			// auto B1 = wB.at(bestPair.second);
			// auto B2 = wB.at(bestPair2.second);
			// wB.removeAll(B1);
			// wB.removeAll(B2);
			continue;
		}
//		qDebug() << "Matching point : " << minDist << "\n" << wA.at(bestPair.first)  << "\n" << wB.at(bestPair.second);
		matches << match;
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

template<typename T>
cimg_library::CImg<T> PhotoAnalysis::fondamentale(QList<QList<T> > matches)
{
	//on construit une matrice A de taille n×9 où chaque ligne est de la forme :
	// [uu′,uv′,u,vu′,vv′,v,u^′,v^′,1]
	CImg<T> A(matches.size(), 9);
	for (int i = 0; i < matches.size();++i) {
		const auto & m=matches.at(i);
		A(i, 0) = m.at(0)*pow(2, m.at(2))*m.at(3)*pow(2, m.at(5));
		A(i, 1) = m.at(0)*pow(2, m.at(2))*m.at(4)*pow(2, m.at(5));
		A(i, 2) = m.at(0)*pow(2, m.at(2));
		A(i, 3) = m.at(1)*pow(2, m.at(2))*m.at(3)*pow(2, m.at(5));
		A(i, 4) = m.at(1)*pow(2, m.at(2))*m.at(4)*pow(2, m.at(5));
		A(i, 5) = m.at(1)*pow(2, m.at(2));
		A(i, 6) = m.at(3)*pow(2, m.at(5));
		A(i, 7) = m.at(4)*pow(2, m.at(5));
		A(i, 8) = 1;
	}
	auto r = A.get_SVD();
	// qDebug() << r.at(0).width() << r.at(0).height();
	// qDebug() << r.at(1).width() << r.at(1).height();
	// qDebug() << r.at(2).width() << r.at(2).height();
	// La solution f est la dernière colonne de V
	// On extrait cette colonne (qui contient n éléments)
	CImg<T> f = r.at(2).get_column(r.at(2).height() - 1);
	// On ne garde que les 9 premiers coefficients (la partie utile)
	// et on les remet en forme 3x3 pour former F
	CImg<T> F(3, 3);
	for (int i = 0; i < 9; ++i) {
		F(i % 3, i / 3) = f(i); // Attention à l'ordre de remplissage
	}
	CImg<T> U, S, V;
	F.SVD(U, S, V);

	// Forcer la plus petite valeur singulière à zéro
	int idx_min = 0;
	for (int i = 0; i < (int)S.size(); ++i)
		if (S(i) < S(idx_min)) idx_min = i;
	S(idx_min) = 0;
	S.diagonal();

	// Reconstruction F_corrigée = U * diag(S) * V^T
	CImg<T> F_corr = U * S * V.get_transpose();
	cimg_forX(F_corr, x) {
		qDebug() << F_corr(x, 0) << F_corr(x, 1) << F_corr(x, 2);
	}
	return F_corr;
}

template<typename T>
CImg<T> PhotoAnalysis::placePoint(const CImg<T> &P1, const CImg<T> &P2, CImg<T> x1, CImg<T> x2) {
	// P1, P2 : CImg<double> 3x4 (les matrices de projection)
	// x1, x2 : CImg<double> 3x1 (points homogènes, ou 1x3)
	CImg<T> A(4, 4, 1, 1, 0);

	auto remplir = [&](const CImg<T>& P, const CImg<T>& x, int ligne) {
		T u = x(0), v = x(1);
		for (int j = 0; j < 4; ++j) {
			A(ligne,   j) = u * P(2, j) - P(0, j);
			A(ligne+1, j) = v * P(2, j) - P(1, j);
		}
	};
	remplir(P1, x1, 0);
	remplir(P2, x2, 2);

	CImg<T> U, S, V;
	A.SVD(U, S, V);
	CImg<T> Xh = V.get_column(3);   // vecteur homogène 4x1

	// Déshomogénéisation
	double w = Xh(3);
	if (std::abs(w) < 1e-12) { /* point à l'infini : direction pure */ }
	CImg<T> X = Xh.get_rows(0, 3) / w;   // coordonnées 3D
	return X;
}

void PhotoAnalysis::export2Ply(QList <QList <float>> points3D) {
	QFile file("heightMap.ply");
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
		return;
	}
	QTextStream out(&file);
	out << "ply" << Qt::endl << "format ascii 1.0" << Qt::endl << "element vertex " << points3D.size() << Qt::endl;
	out << "property float x" << Qt::endl << "property float y" << Qt::endl << "property float z" << Qt::endl << "end_header" << Qt::endl;
	for (const auto &p : points3D) {
		out << p.at(0) << "\t" << p.at(1) << "\t" << p.at(2) << Qt::endl;
	}
}