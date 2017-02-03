#include "GraphImplementation.h"

void GraphImplementation::addEdge(UINT n0, UINT n1)
{
	base_graph::addEdge(n0, n1);
}

void GraphImplementation::addTri(UINT n0, UINT n1, UINT n2)
{
	base_graph::addTri({ n0,n1,n2 });
}

void GraphImplementation::addSqr(UINT n0, UINT n1, UINT n2, UINT n3)
{
	base_graph::addSq({ n0, n1, n2, n3 });
}

void GraphImplementation::addTet(UINT n0, UINT n1, UINT n2, UINT n3)
{
	base_graph::addTet({ n0, n1, n2, n3 });
}

void GraphImplementation::addPyr(UINT n0, UINT n1, UINT n2, UINT n3, UINT n4)
{
	base_graph::addPyr({ n0, n1, n2, n3, n4 });
}

void GraphImplementation::addWedge(UINT n0, UINT n1, UINT n2, UINT n3, UINT n4, UINT n5)
{
	base_graph::addWedge({ n0, n1, n2, n3, n4, n5 });
}

void GraphImplementation::addHexa(UINT n0, UINT n1, UINT n2, UINT n3, UINT n4, UINT n5, UINT n6, UINT n7)
{
	base_graph::addHexa({ n0, n1, n2, n3, n4, n5, n6, n7 });
}
