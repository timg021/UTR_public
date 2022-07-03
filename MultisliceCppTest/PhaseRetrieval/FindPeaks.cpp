#include <omp.h>
#include <vector>

#include "XA_ini.h"
#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray3D.h"
#include "XA_data.h"

using namespace xar;

int FindPeaks(XArray3D<float> xar3D, double datomsizeXY, double datomsizeZ, int natommax, string XYZfilename)
// finds peaks (maximums) in the distribution xar3D, with no more than 1 peak per parallelepiped with dimensions (datomsizeXY x datomsizeZ)
// natommax - max number of peaks to be saved
// XYZfilename - name of a file for saving the found peaks in the Kirkland's XYZ format
{
	const Wavehead3D* phead3D = dynamic_cast<const Wavehead3D*>(xar3D.GetHeadPtr());
	if (phead3D == nullptr) throw std::runtime_error("Error: input 3D array does not have a header in FindPeaks().");

	int nx = (int)xar3D.GetDim3();
	int ny = (int)xar3D.GetDim2();
	int nz = (int)xar3D.GetDim1();
	double xlo = phead3D->GetXlo();
	double ylo = phead3D->GetYlo();
	double zlo = phead3D->GetZlo();
	double xhi = phead3D->GetXhi();
	double yhi = phead3D->GetYhi();
	double zhi = phead3D->GetZhi();
	double xst = phead3D->GetXStep(nx);
	double yst = phead3D->GetXStep(ny);
	double zst = phead3D->GetXStep(nz);

	// finding atomic positions

	// exclude the points located outside the reconstruction volume from the subsequent peak search
	double xR = (xhi - xlo) / 2.0; // x-radius
	double yR = (yhi - ylo) / 2.0; // y-radius
	double zR = (zhi - zlo) / 2.0; // z-radius
	double R2 = (xR * xR + yR * yR + zR * zR) / 3.0;
	float K3min = (float)xar3D.Norm(eNormMin);

	#pragma omp parallel for shared(xar3D, K3min, xR, yR, zR, R2, xst, yst, zst)
	for (int k = 0; k < nz; k++)
	{
		double zzz = -zR + zst * k;
		zzz *= zzz;
		for (int j = 0; j < ny; j++)
		{
			double yyy = -yR + yst * j;
			yyy *= yyy; yyy += zzz;
			for (int i = 0; i < nx; i++)
			{
				double xxx = -xR + xst * i;
				xxx *= xxx; xxx += yyy;
				if (xxx > R2)
					xar3D[k][j][i] = K3min; // mark points outside the reconstruction volume for exclusion from the peak search
			}
		}
	}

	int katom = int(datomsizeZ / zst + 0.5), jatom = int(datomsizeXY / yst + 0.5), iatom = int(datomsizeXY / xst + 0.5);
	int natom(0);
	vector<int> vimax, vjmax, vkmax;
	vector<float> xa, ya, za;

	// search for peaks
	// NOTE that we exclude one-atomsize-wide vicinity of the outer boundary from the search, as we expect artefacts there
	#pragma omp parallel for shared(xar3D, natom, katom, jatom, iatom, vimax, vjmax, vkmax)
	for (int k = katom; k < nz - katom * 2; k += katom)
	{
		for (int j = jatom; j < ny - jatom * 2; j += jatom)
		{
			for (int i = iatom; i < nx - iatom * 2; i += iatom)
			{
				double K3max(0);
				int kmax(0), jmax(0), imax(0);
				for (int kk = k; kk < k + katom; kk++)
					for (int jj = j; jj < j + jatom; jj++)
						for (int ii = i; ii < i + iatom; ii++)
							if (xar3D[kk][jj][ii] > K3max)
							{
								K3max = xar3D[kk][jj][ii];
								xar3D[kmax][jmax][imax] = 0.0;
								kmax = kk; jmax = jj; imax = ii;
							}
							else xar3D[kk][jj][ii] = 0.0;
				#pragma omp critical
				if (K3max > 0)
				{
					natom++;
					vimax.push_back(imax);
					vjmax.push_back(jmax);
					vkmax.push_back(kmax);
				}
			}
		}
	}

	if (natom == 0) printf("\n\n!!WARNING: no peaks have been found!");
	else
	{

		// create vectors of peak location coordinates
		xa.resize(natom); ya.resize(natom); za.resize(natom);
		vector<Pair2> K3maxPair(natom);
		for (int n = 0; n < natom; n++)
		{
			xa[n] = float(xlo + xst * vimax[n]);
			ya[n] = float(ylo + yst * vjmax[n]);
			za[n] = float(zst * vkmax[n]);
			K3maxPair[n].v = xar3D[vkmax[n]][vjmax[n]][vimax[n]];
			K3maxPair[n].n = n;
		}

		// sort the located peaks according to their height, in order to improve the subsequent procedure of elimination of closely located peaks
		printf("\nSorting the located peaks in descending order according to their heights ...");
		std::qsort(&K3maxPair[0], (size_t)natom, sizeof(Pair2), Pair2comp);

		// exclude the smaller one from each pair of peak positions that are located closer than datomsize to each other (e.g. in adjacent corners of neigbouring cubes)
		double datomsize2 = datomsizeXY * datomsizeXY;
		vector<int> vimax1, vjmax1, vkmax1, Znum1;
		vector<float> xa1, ya1, za1, occ1, wobble1;
		printf("\nEliminating adjacent peaks in the reconstructed 3D distribution ...");
		int n, m;
		for (int nn = natom - 1; nn >= 0; nn--)
		{
			if (K3maxPair[nn].v != 0.0)
			{
				// we can already count this maximum in, as it is guaranteed to be larger than all subsequent ones
				n = K3maxPair[nn].n;
				Znum1.push_back(6);
				vimax1.push_back(vimax[n]);
				vjmax1.push_back(vjmax[n]);
				vkmax1.push_back(vkmax[n]);
				xa1.push_back(xa[n]);
				ya1.push_back(ya[n]);
				za1.push_back(za[n]);
				occ1.push_back((float)K3maxPair[nn].v);
				wobble1.push_back(0);

				#pragma omp parallel for private(m)
				for (int mm = nn - 1; mm >= 0; mm--)
				{
					if (K3maxPair[mm].v != 0.0)
					{
						m = K3maxPair[mm].n;
						if ((xa[n] - xa[m]) * (xa[n] - xa[m]) + (ya[n] - ya[m]) * (ya[n] - ya[m]) + (za[n] - za[m]) * (za[n] - za[m]) < datomsize2)
							K3maxPair[mm].v = 0.0;
					}
				}
			}
		}
		natom = (int)Znum1.size();
		printf("\n%d peak positions have been found in total.", natom);

		// copy the located isolated peaks back into the input array
		xar3D.Fill(0.0);
		for (int n = 0; n < natom; n++)
			xar3D[vkmax1[n]][vjmax1[n]][vimax1[n]] = occ1[n];

		if (natom > natommax)
			printf("\nNot saving the reconstructed peak-localized 3D object into an XYZ file, since too many peaks (%d) have been found.", natom);
		else
		{
			size_t dotpos = XYZfilename.find_last_of(".");
			string fileoutXYZ = XYZfilename.replace(dotpos + 1, XYZfilename.size() - dotpos - 1, "xyz");
			printf("\nSaving the reconstructed peak-localized 3D object into an XYZ file %s ...", fileoutXYZ.c_str());
			SaveXYZfile(fileoutXYZ, float(xhi - xlo), float(zhi - zlo), (int)natom, &Znum1[0], &xa1[0], &ya1[0], &za1[0], &occ1[0], &wobble1[0]);
		}
	}

	return natom;
}