/*
 *  This file is part of birch-clustering-algorithm.
 *
 *  birch-clustering-algorithm is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  birch-clustering-algorithm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with birch-clustering-algorithm.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Copyright (C) 2011 Taesik Yoon (otterrrr@gmail.com)
 */

#ifndef __CFTREE_REDIST_H__
#define __CFTREE_REDIST_H__

// class CFTree
// {

	private:
		bool _has_differences( std::vector<ublas_vec_type>& lhs, std::vector<ublas_vec_type>& rhs)
		{
			assert( lhs.size() == rhs.size() );

			for( std::size_t i = 0 ; i < lhs.size() ; i++ )
			{
				if( norm_2(lhs[i] - rhs[i]) > std::numeric_limits<float_type>::epsilon() )
					return true;
			}
			return false;
		}

		bool _has_differences( std::vector<ublas_vec_type>& lhs, std::vector<ublas_vec_type>& rhs, std::vector<bool>& active )
		{
			assert( lhs.size() == rhs.size() );

			for( std::size_t i = 0 ; i < lhs.size() ; i++ )
				active[i] = norm_2(lhs[i] - rhs[i]) > std::numeric_limits<float_type>::epsilon();

			return std::count( active.begin(), active.end(), true ) > 0;
		}

	public:

		template<typename item_list_type>
		void redist_kmeans( item_list_type& items, cfentry_vec_type& entries, std::size_t iteration = 2 )
		{
			using namespace boost::numeric::ublas;
			
			if( items.empty() )
				return;

			assert(items[0].size() == dim);

			// start from k means from k entries
			std::vector<ublas_vec_type> prev_means(entries.size());
			std::vector<ublas_vec_type> means(entries.size());
			for( std::size_t i = 0 ; i < means.size() ; i++ )
			{
				prev_means[i].resize( dim );
				prev_means[i].clear();

				CFEntry& e = entries[i];
				ublas_vec_type& mean = means[i];

				mean.resize( dim );
				std::copy( e.sum, e.sum + dim, mean.begin() );
				mean /= e.n;
			}

			// until it is converged
			std::size_t iteration_count = 0;
			if( iteration == 0 )
				iteration = (std::numeric_limits<std::size_t>::max)();

			bool active = true;

			//while( iteration_count < iteration && _has_differences( prev_means, means ) )
			while( iteration_count < iteration && active )
			{
				active = false;
				for( item_list_type::iterator item_it = items.begin() ; item_it != items.end() ; ++item_it )
				{
					item_list_type::value_type& item = *item_it;
					float_type min_dist = (std::numeric_limits<float_type>::max)();

					int prev_cid = item.cid();

					for( std::size_t cid = 0 ; cid < means.size() ; ++cid )
					{
						ublas_vec_type diff(dim);
						std::transform( &item[0], &item[0] + dim, means[cid].begin(), diff.begin(), std::minus<float_type>() );
						float_type dist = norm_2( diff );

						if( min_dist > dist )
						{
							min_dist = dist;
							item_it->cid() = cid;
						}
					}

					if( prev_cid != item.cid() )
						active = true;
				}

				std::stringstream ss;
				ss << "k-means_iteration" << iteration_count << ".txt";
				std::ofstream fout( ss.str().c_str() );

				for( std::size_t c = 0 ; c < means.size() ; c++ )
				{
					fout << "(" << c << ") ";
					for( std::size_t d = 0 ; d < dim ; d++ )
						fout << means[c][d] << (d == dim-1 ? "" : "," );
					fout << std::endl;
				}
				fout << std::endl;

				for( std::size_t i = 0 ; i < items.size() ; i++ )
					fout << i << ":" << items[i].cid() << std::endl;
				fout.close();

				// store means to prev_means and zeoring means
				prev_means = means;
				for( std::size_t i = 0 ; i < means.size() ; i++ )
					means[i].clear();

				// rearrange means and count # items for each cluster
				std::vector<std::size_t> mean_counts(means.size(), 0);
				for( item_list_type::iterator item_it = items.begin() ; item_it != items.end() ; ++item_it )
				{
					item_list_type::value_type& item = *item_it;
					std::transform( &item[0], &item[0] + dim, means[item.cid()].begin(), means[item.cid()].begin(), std::plus<float_type>() );
					++mean_counts[item.cid()];
				}

				// averaging means to calculate centroids
				for( std::size_t i = 0 ; i < means.size() ; i++ )
					means[i] /= mean_counts[i];

				iteration_count++;
			}
			//std::cout << "iteration count = " << iteration_count << std::endl;
		}

		/************************************************************************/
		/* The original redistribution code of birch
		/* In my view point, it could be burdensome due to O(n^2) cost
		/************************************************************************/
		struct subcluster_summary
		{
			subcluster_summary( ) : radius(0.0), norm(0.0) {}
			subcluster_summary( const ublas_vec_type& in_center, const float_type& in_radius, const float_type& in_norm ) : center( in_center ), radius(in_radius), norm(in_norm) {}
			//subcluster_summary( float_type* in_center, float_type& in_radius, float_type& in_norm ) : center( in_center, in_center + dim ), radius(in_radius), norm(in_norm) {}

			ublas_vec_type	center;
			float_type		radius;
			float_type		norm;
		};

		typedef std::vector< subcluster_summary > subsum_vec_type;

		struct subcluster_lessthan_norm
		{
			bool operator() ( const subcluster_summary& lhs, const subcluster_summary& rhs ) const { return (lhs.norm) < (rhs.norm); }
		};

		template<typename _iter>
		void redist( _iter begin, _iter end, cfentry_vec_type& entries, std::vector<int>& out_cid )
		{
			using namespace boost::numeric::ublas;

			// prepare summaries for each subcluster
			// summaries = ( center, radius, norm )
			subsum_vec_type subclusters;
			subclusters.reserve(entries.size());
			for( std::size_t i = 0 ; i < entries.size() ; i++ )
			{
				const CFEntry& e = entries[i];
				ublas_vec_type center(dim);
				std::copy(e.sum, e.sum + dim, center.begin());
				center /= e.n;
				subclusters.push_back( subcluster_summary( center, _Radius(e), std::sqrt(inner_prod(center, center) )) ); 
			}

			std::sort( subclusters.begin(), subclusters.end(), subcluster_lessthan_norm() );

			// in addition to an individual summary for each subcluster
			// calculate pairwise euclidean distances of subclusters
			std::size_t n = subclusters.size();
			ublas_sym_matrix_type dist_mat(n,n);
			for( std::size_t i = 0 ; i < n-1 ; i++ )
			{
				for( std::size_t j = i+1 ; j < n ; j++ )
				{
					ublas_vec_type diff = subclusters[i].center - subclusters[j].center;
					dist_mat(i,j) = inner_prod(diff, diff);
				}
			}

			out_cid.clear();
			out_cid.reserve(end - begin);
			for( _iter it = begin ; it != end ; it++ )
			{
				ublas_vec_type v(dim);
				std::copy(&(*it)[0], &(*it)[0] + dim, v.begin());
				out_cid.push_back( _redist( v, subclusters, dist_mat ) );
			}
		}

	private:

		int _redist( ublas_vec_type& tmpv, subsum_vec_type& subsums, ublas_sym_matrix_type& dist_mat )
		{
			int    imin,imax,i,k,n,start,end,median;
			float d,tmpnorm,idist,tmpdist;
			ublas_vec_type diff;

			i = 0;
			n = (int)subsums.size() ;
			tmpnorm = std::sqrt( inner_prod(tmpv, tmpv) );

			// i=ClosestNorm(tmpnorm,norms,0,n-1);
			// for efficiency, replace recursion by iteration
			start=0;
			end=n-1;
			while(start<end)
			{
				if (end-start==1)
				{
					float_type norm_end = subsums[end].norm;
					float_type norm_start = subsums[start].norm;

					i = tmpnorm > norm_end ? end :
						tmpnorm < norm_start ? start :
						tmpnorm - norm_start < norm_end - tmpnorm ? start : end;
					start = end = i;
				}
				else
				{
					median = (start+end)/2;
					float_type norm_med = subsums[median].norm;
					if (tmpnorm > norm_med)
						start=median;
					else
						end=median;
				}
			}

			diff = tmpv - subsums[i].center;
			idist= inner_prod(diff, diff);

			// imin=MinLargerThan(tmpnorm-sqrt(idist),norms,0,n-1);
			// imax=MaxSmallerThan(tmpnorm+sqrt(idist),norms,0,n-1);
			// for efficiency, replace recursion by iteration

			tmpdist=tmpnorm-sqrt(idist);
			start=0;
			end=n-1;
			while (start<end)
			{
				median=(start+end)/2;
				float_type norm_med = subsums[median].norm;
				if (tmpdist > norm_med)
					start=median+1;
				else
					end=median;
			}
			imin=start;

			tmpdist=tmpnorm+sqrt(idist);
			start=0;
			end=n-1;
			while(start<end)
			{
				median=(start+end+1)/2;
				float_type norm_med = subsums[median].norm;
				if (tmpdist < norm_med)
					end=median-1;
				else
					start=median;
			}
			imax=start;

			// ClosestCenter(i,idist,tmpv,centers,imin,imax,matrix,n);
			// for efficiency, replace procedure by inline
			k=imin;
			while (k<=imax)
			{
				if (dist_mat(k,i) <= 4*idist)
				{
					diff = tmpv - subsums[k].center;
					d = inner_prod(diff,diff);
					if (d < idist)
					{
						idist=d;
						i=k;
					}
				}
				k++;
			}
			return i;
		}
// }

#endif
