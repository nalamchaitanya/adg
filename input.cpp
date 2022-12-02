#include<iostream>
#include<vector>
using namespace std;
int main()
{
	int n,m;
	cin >> n >> m;
	int *a = new int[n + 1];
	int *lst = new int[2*m];
	a[0] = n;
	vector<int> g[n + 1];
	for(int i = 0; i < m; i ++)
	{
		int u,v;
		cin >> u >> v;
		g[u].push_back(v);
		g[v].push_back(u);
	}
	int ctr = 0;
	for(int i = 1; i <=n; i ++)
	{
		a[i] = ctr;
		for(auto &x: g[i])
		{
			lst[ctr++] = x;
		}
	}
	for(int i = 0; i <=n ; i++)
	{
		cout << a[i] <<" ";
	}
	cout << endl;
	for(int i = 0; i < m; i ++)
	{
		cout << lst[i] <<" ";
	}
	cout << endl;

}