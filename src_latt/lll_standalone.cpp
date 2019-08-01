#include <vector>
#include <algorithm>
#include <cmath>

typedef int Tint;
typedef float Tfloat;

inline int ind( int i, int j, int n) {
	return i*n+j;
}

// return Greduced, U s.t. Gred = Ut*G*U
std::pair<std::vector<Tint>, std::vector<Tint>> lll( std::vector<Tint> G, int n) {
	Tfloat delta = 0.99;
	Tfloat eta = 0.501;
	int kappa = 1;
	std::vector<Tfloat> mu(n*n, 0);
	std::vector<Tfloat> r(n,0);
	std::vector<Tint> U(n*n,0);

	// init mu, r, U
	for( int i = 0; i < n; i++ ) {
		mu[ind(i,i,n)] = 1.;
		U[ind(i,i,n)] = 1;
	}
	r[0] = G[ind(0,0,n)];
	int rnd = 0;
	while( kappa < n ) {

		// compute mu
		for( int j = 0; j < kappa; j++) {
			mu[ind(kappa, j, n)] = G[ind(kappa, j, n)];
			for( int i = 0; i < j; i++ )
				mu[ind(kappa,j,n)] -= mu[ind(kappa,i,n)] * r[i] * mu[ind(j,i,n)];
			mu[ind(kappa,j,n)] /= r[j];
		}

		// compute r[kappa]
		r[kappa] = G[ind(kappa, kappa, n)];
		for( int i = 0; i < kappa; i++)
			r[kappa] -= mu[ind(kappa,i,n)] * r[i] * mu[ind(kappa, i,n)];

		// eta-size reduce
		for( int i = kappa-1; i >= 0; i--) {
			if( std::fabs(mu[ind(kappa,i,n)]) > eta) {
				rnd = std::lround(mu[ind(kappa,i,n)]);
			
				// update G
				G[ind(kappa, kappa, n)] -= 2*rnd*G[ind(kappa,i,n)] - rnd*rnd*G[ind(i,i,n)];
				for( int j = 0; j < kappa; j++ ) {
					G[ind(kappa,j,n)] -= rnd * G[ind(i,j,n)];
					G[ind(j,kappa,n)] = G[ind(kappa,j,n)];
				}
				for( int j = kappa+1; j < n; j++) {
					G[ind(kappa,j,n)] -= rnd * G[ind(i,j,n)];
					G[ind(j,kappa,n)] = G[ind(kappa,j,n)];
				}

				// update U
				for( int j = 0; j < n; j++ ) {
					U[ind(j,kappa,n)] -= rnd * U[ind(j,i,n)];
				}

				for( int j = 0; j <= i; j++ )
					mu[ind(kappa,j,n)] -= rnd * mu[ind(i,j,n)];
			}
		}

		// lovasz condition
		if( (delta - mu[ind(kappa,kappa-1, n)]*mu[ind(kappa, kappa-1, n)]) * r[kappa-1] <= r[kappa] ) {
			kappa++;
		}
		else {
			// swap kappa-1 and kappa
			for( int j = 0; j < n; j++ ) 
				std::swap( G[ind(kappa-1,j,n)], G[ind(kappa,j,n)] );
			for( int i = 0; i < n; i++ ) 
				std::swap( G[ind(i,kappa-1,n)], G[ind(i,kappa,n)] );
			if( kappa == 1 )
				r[0] = G[ind(0,0,n)];
			// update U
			for( int i = 0; i < n; i++ ) 
				std::swap( U[ind(i,kappa-1, n)], U[ind(i, kappa, n)]);

			kappa = std::max(1, kappa-1);
		}
	}
	return {G, U};
}