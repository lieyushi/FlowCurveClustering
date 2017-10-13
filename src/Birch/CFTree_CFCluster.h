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

#ifndef	__CFCLUSTER_H__
#define __CFCLUSTER_H__


/************************************************************************/
/* a partial class of CFTree, CFCluster
/************************************************************************/

// class CFTree
// {

	public:
		void cluster( cfentry_vec_type& entries )
		{
			get_entries(entries);
			_cluster(entries);
		}

	private:

		struct HierarchicalClustering
		{
			typedef boost::numeric::ublas::symmetric_matrix<float_type> dist_matrix_type;

			HierarchicalClustering(int n, dist_func_type& in_dist_func) : size(n), step(-1), ii(n), jj(n), cf(n), dd(n), dist_func(in_dist_func), chain(n+1), chainptr(-1), stopchain(FALSE)
			{
			}

			void merge( cfentry_vec_type& entries )
			{
				int		nentry = (int)entries.size();
				int 	i,j,n1,n2;

				int		CurI, PrevI, NextI;
				int 	uncheckcnt = nentry;
				
				std::vector<int> checked(nentry);
				for (i=0;i<nentry;i++)
					checked[i]=i+1;
				// 0: invalid after being merged to other entries
				// positive 1..nentry+1 :     original entries
				// negative -1..-(nentry-1) : merged entries

				// get initial distances
				// std::vector<float_type> dist(nentry*(nentry-1)/2);
				dist_matrix_type dist(nentry, nentry);

				for (i=0; i<nentry-1; i++)
					for (j=i+1; j<nentry; j++)
						dist(i, j) = dist_func(entries[i],entries[j]);

				CurI = rand() % nentry;			// step1 
				chain[++chainptr]=CurI;

				while (uncheckcnt>1)
				{
					// step4
					if (chainptr==-1)
					{
						chainptr++; 
						chain[chainptr]=pick_one_unchecked(nentry, &checked[0]);
					}
					PrevI = chainptr > 0 ? chain[chainptr-1] : -1;
					stopchain=FALSE;

					// step2
					while (stopchain==FALSE)
					{
						CurI=chain[chainptr];
						NextI = nearest_neighbor(CurI,nentry,&checked[0], dist);
						
						// it is impossible NextI be -1 because uncheckcnt>1
						if (NextI==PrevI)
							stopchain = TRUE;
						else
						{
							chain[++chainptr]=NextI;
							PrevI = CurI;
						}
					} // end of while for step 2

					step++;
					
					// step3
					ii[step] = checked[CurI];
					jj[step] = checked[NextI];

					dd[step] = dist(CurI, NextI);

					bool curr_org = checked[CurI] > 0;
					bool next_org = checked[NextI] > 0;
					
					CFEntry& curr_entry = curr_org ? entries[CurI] : cf[-checked[CurI]-1];
					CFEntry& next_entry = next_org ? entries[NextI] : cf[-checked[NextI]-1];
					n1 = curr_entry.n;
					n2 = next_entry.n;
					cf[step] = curr_entry + next_entry;

					update_distance(n1,n2,CurI,NextI,nentry,&checked[0], dist);
					uncheckcnt--;
					checked[CurI] = -(step+1);	    
					checked[NextI] = 0;
					chainptr--;
					chainptr--;
				} //end of while (uncheckcnt>1)

				// prepare for SplitHierarchy
				stopchain = FALSE;
				chainptr = 0;
				chain[chainptr] = -(step+1);
			}

			void split( float_type ft )
			{
				while( _split(ft) )
					/* nothing to do */;
			}

			short _split( float_type ft )
			{
				int i, j;
				
				if ( chainptr == size )
					return FALSE;

				i = largest_merge(chainptr, ft);
				if (i!=-1)
				{
					j = -chain[i]-1;
					chain[i] = ii[j];
					chain[++chainptr]=jj[j];
					return TRUE;
				}
				
				stopchain = TRUE;
				return FALSE;
			}

			void result( /* inout */cfentry_vec_type& entries )
			{
				int j;
				std::vector<CFEntry> tmpentries( chainptr + 1 );
				for( j = 0 ; j <= chainptr ; j++ )
				{
					tmpentries[j] = chain[j] < 0 ? cf[-chain[j]-1] : entries[chain[j]-1];
				}
				entries = tmpentries;
			}

			/* for SplitHierarchy use only */
			int farthest_merge(int chainptr) 
			{
				if (chainptr<=0)
					return chainptr;

				float d, dmax = 0;
				int    i, imax = -1;
				for (i=0; i<=chainptr; i++)
				{
					if (chain[i]<0)
					{
						d = dd[-chain[i]-1];
						if (d>dmax) {imax = i; dmax = d;}
					}
				}
				return imax;
			}

			/* for SplitHierarchy use only */
			int largest_merge(int chainptr, float_type dist_threshold)
			{
				for (int i=0; i<=chainptr; i++)
				{
					if (chain[i]<0)
					{
						if ( _Diameter(cf[-chain[i]-1]) > dist_threshold)
							return i;
					}
				}
				return -1;
			}

			/* for MergeHierarchy use only */
			int nearest_neighbor(int CurI, int n, int *checked, dist_matrix_type& dist)
			{
				int    imin=0;
				float d, dmin = (std::numeric_limits<float_type>::max)();
				for( int i = 0 ; i < n ; i++ )
				{
					if( i == CurI || checked[i] == 0 )
						continue;

					d = dist(i, CurI);
					if (d < dmin)
					{
						dmin=d;
						imin=i;
					}
				}

				return dmin < (std::numeric_limits<float_type>::max)() ? imin : -1;
			}

			/* for MergeHierarchy use only */
			int pick_one_unchecked(int n, int *checked)
			{
				int i,j = rand() % n;
				for (i=0;i<n;i++) 
					if (checked[(i+j)%n]!=0)
						return (i+j)%n;
				return -1;
			}

			/* for MergeHierarchy use only */
			void update_distance(int n1, int n2, int CurI, int NextI, int n, int *checked, dist_matrix_type& dist)
			{
				for( int i = 0 ; i < n ; i++ )
				{
					if( i == CurI || i == NextI )
						continue;

					if( checked[i] != 0 )
						dist(i, CurI) = (n1 * dist(i, CurI) + n2 * dist(i, NextI)) / (n1 + n2);
				}
			}

			int						size;
			int						step;
			std::vector<int>		ii;
			std::vector<int>		jj;
			std::vector<CFEntry>	cf;
			std::vector<float_type>	dd;
			
			std::vector<int>		chain;
			int						chainptr;
			short					stopchain;
			dist_func_type&			dist_func;
		};

		void refine_cluster( cfentry_vec_type& entries )
		{
			std::vector<bool> merged(entries.size(), false);

			std::vector<std::size_t> not_visited;
			not_visited.reserve( entries.size() );

			for( std::size_t i = 0 ; i < entries.size() ; i++ )
				not_visited.push_back(i);
			
			while( !not_visited.empty() )
			{
				// pick any entry
				CFEntry& ref_entry = entries[not_visited.back()];
				CFEntry curr_entry = ref_entry;
				not_visited.pop_back();

				bool something_merged = false;
				for( std::size_t i = 0 ; i < not_visited.size() ; i++ )
				{
					// index of next item
					std::size_t v = not_visited[i];
					CFEntry& e = entries[ v ];
					if( dist_func(ref_entry, e) <= dist_threshold/2 )
					{
						curr_entry += e;
						merged[ v ] = true;
						something_merged = true;
					}
				}

				if( something_merged )
					ref_entry = curr_entry;

				// remove if visited
				struct _remove_if_merged
				{
					_remove_if_merged( std::vector<bool>& in_visited ) : visited(in_visited) {}
					bool operator()( const std::size_t i )  const { return visited[i]; }
				private:
					std::vector<bool>& visited;
				};

				// prepare for next loop, removing visited indices
				if( something_merged )
					not_visited.erase( std::remove_if( not_visited.begin(), not_visited.end(), _remove_if_merged(merged)), not_visited.end() );
			}

			struct _remove_if_merged_by_item
			{
				typedef CFEntry argument_type;
				_remove_if_merged_by_item( cfentry_vec_type& in_entries, std::vector<bool>& in_merged ) : entries(in_entries), merged(in_merged) {}
				bool operator() ( const CFEntry& e ) const
				{
					std::ptrdiff_t i = (&e - &entries[0]);
					return merged[i];
				}
				
			private:
				std::vector<bool>& merged;
				cfentry_vec_type& entries;
			};
			entries.erase( std::remove_if( entries.begin(), entries.end(), _remove_if_merged_by_item( entries ,merged) ), entries.end());
		}

		void _cluster( cfentry_vec_type& entries )
		{
			int n = (int)entries.size();

			if( n <= 1 )
				return;

			if( dist_func == _DistD0 || dist_func == _DistD1 )
			{
				refine_cluster( entries );
			}
			else
			{
				HierarchicalClustering h( n - 1, dist_func );
				h.merge( entries );
				h.split( dist_threshold );
				h.result( entries );
			}
		}

// };


#endif
