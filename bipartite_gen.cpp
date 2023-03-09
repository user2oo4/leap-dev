/*
khoi orz, go check out his algo
-normie-
*/
#include <bits/stdc++.h>
using namespace std;
#define rep(i,n) for(int64_t i=0;i < (int64_t)(n);i++)
#pragma comment(linker, "/stack:200000000")
#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#define FILE_IN "ninja.inp"
#define FILE_OUT "ninja.out"
#define ofile freopen(FILE_IN,"r",stdin);freopen(FILE_OUT,"w",stdout)
#define fio ios::sync_with_stdio(0);cin.tie(0);cout.tie(0)
#define nfio cin.tie(0);cout.tie(0)
#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))
#define ord(a,b,c) ((a>=b)and(b>=c))
#define MOD (ll(1000000007))
#define MAX 300001
#define mag 1048576
#define fi first
#define se second
#define pow2(x) (ll(1)<<x)
#define pii pair<int,int>
#define piii pair<int,pii>
#define pll pair<ll,ll>
#define plll pair<ll,pll>
#define For(i,__,___) for(int i=__;i<=___;i++)
#define Rep(i,__,___) for(int i=__;i>=___;i--)
#define ordered_set tree<long long,null_type,less<long long>,rb_tree_tag,tree_order_statistics_node_update>
#define endl "\n"
#define bi BigInt
#define ll long long
#define pi 3.1415926535897
//------START-----------//
struct dsu
{
	int n,par[300100],h[300100],sz[300100];
	dsu (int n=0)
	{
		for (int i=1;i<=n;i++)
		{
			par[i]=i;
			sz[i]=1;
			h[i]=1;
		}
		this->n=n;
	}
	
	void reset (int n)
	{
		for (int i=1;i<=n;i++)
		{
			par[i]=i;
			sz[i]=1;
			h[i]=1;
		}
		this->n=n;
	}
	int get_par(int x)
	{
		if (par[par[x]]==par[x]) return par[x];
		else return par[x]=get_par(par[x]);
	}
	int check_same(int a, int b)
	{
		return (get_par(a)==get_par(b));
	}
	int add_edge(int a, int b)
	{
		int ha=get_par(a),hb=get_par(b);
		if (ha!=hb)
		{
			if (h[ha]<h[hb])
			{
				par[ha]=hb;
				sz[hb]+=sz[ha];
				return hb;
			}
			else
			if (h[ha]>h[hb])
			{
				par[hb]=ha;
				sz[ha]+=sz[hb];
				return ha;
			}
			else
			{
				par[hb]=ha;
				sz[ha]+=sz[hb];
				h[ha]++;
				return ha;
			}
			return 1;
		}
		else return 0;
	}
};
//------END-----------//
ll n,m,k,t,t1,i,j,a,b,u,res;
ll val[101];
int main()
{
    srand(95746);
    fio;
    
	cin>>n;
    cout<<n<<endl;
	
	for (i=0;i<n;i++) {
        u=rand()%101;
		u+=50;
		
        u*=(rand()&1)*2-1;
		cout<<u<<endl;
	}

    for (i=0;i<n/2;i++) for (j=n/2;j<n;j++) {
        u=rand()%101;
		u+=50;
		
        u*=(rand()&1)*2-1;
		
		cout<<i<<' '<<j<<' '<<u<<endl; 
    }
    cout<<-1<<endl;
    // cout<<res<<endl;
}
// a;