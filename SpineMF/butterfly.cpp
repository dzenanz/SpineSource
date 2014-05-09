#include "declarations.h"
#include <Eigen/Eigen>

#ifdef _MSC_VER
inline bool isnan(double x) {
    return x != x;
}
#endif

void topologicallySubdivide(Cell *qe, unsigned level)
{
	CellVertexIterator cvi(qe);
	Vertex *v;
	while ((v=cvi.next())!=0)
	{
		VertexEdgeIterator vei(v);
		Edge *e;
		while ((e=vei.next())!=0) //calculate midpoints
			if (isnan(e->midpoint[0]))
			{
				e->midpoint=(e->Org()->pos + e->Dest()->pos)*0.5;
				e->Sym()->midpoint=e->midpoint;
			}
			else
				; //midpoint has already been calculated by iterating over a different vertex
	}
	
	//now do actual topological face splitting
	CellFaceIterator cfi(qe);
	Face *f;
	//unsigned i=0;
	while ((f=cfi.next())!=0)
	{
		Vertex *a, *b, *c;
		Edge *e=f->getEdge();
		Edge *orig=e;
		
		//get a
		if (isnan(e->midpoint[0]))//this is a new edge
		{
			if (e->Org()->level==level)
				a=e->Org();
			else
			{
				a=e->Dest();
				e=e->Lnext();
			}
		}
		else
		{
			orig=split(e);
			a=e->Org();
			a->level=level;
			a->pos=e->midpoint;
			a->valence=6;
			e->midpoint=vec3(std::numeric_limits<float>::quiet_NaN(),0,0);//deinitialize
			e->Sym()->midpoint=e->midpoint;
		}
		e=e->Lnext();
		//get b
		if (isnan(e->midpoint[0]))//this is a new edge
		{
			b=e->Dest();
			e=e->Lnext();
		}
		else
		{
			split(e);
			b=e->Org();
			b->level=level;
			b->pos=e->midpoint;
			b->valence=6;
			e->midpoint=vec3(std::numeric_limits<float>::quiet_NaN(),0,0);//deinitialize
			e->Sym()->midpoint=e->midpoint;
		}
		e=e->Lnext();
		//get c
		if (isnan(e->midpoint[0]))//this is a new edge
			c=e->Dest();
		else
		{
			split(e);
			c=e->Org();
			c->level=level;
			c->pos=e->midpoint;
			c->valence=6;
			e->midpoint=vec3(std::numeric_limits<float>::quiet_NaN(),0,0);//deinitialize
			e->Sym()->midpoint=e->midpoint;
		}

		//now split face
		qe->makeFaceEdge(f, c, a);
		qe->makeFaceEdge(f, a, b);
		qe->makeFaceEdge(f, b, c);
		//objWriteCell(qe, ("D:\\Temp\\subdivF"+QString::number(++i).rightJustified(2,'0')+".obj").toStdString().c_str());
	}
	//objWriteCell(qe, "D:\\butterfly.obj");
}

Edge* stepAlong(Edge *e, unsigned steps) //in regular mesh
{
	for (unsigned i=1; i<steps; i++)
	{
		e=e->Sym();
		e=e->Onext()->Onext()->Onext();
	}
	return e;
}

void calcSubdivPosExtra3(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	vec3 pos=e->Org()->subdivPos*0.75; //+9/12
	pos+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.4166666666666666666; //+5/12
	e=e->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0833333333333333333; //-1/12
	e=e->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0833333333333333333; //-1/12
	v->subdivPos=pos;
}

void calcSubdivPosExtra4(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	vec3 pos=e->Org()->subdivPos*0.75; //+6/8
	pos+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.375; //+3/8
	e=e->Onext()->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.125; //-1/8
	v->subdivPos=pos;
}

void calcSubdivPosExtraRest(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	vec3 pos=e->Org()->subdivPos*0.75;
	double tc;
	unsigned valence=left->Org()->valence;
	for (unsigned i=0; i<valence; i++)
	{
		tc=1.0/valence*(0.25+cos(2*i*M_PI/valence)+0.5*cos(4*i*M_PI/valence));
		pos+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*tc;
		e=e->Onext();
	}
	v->subdivPos=pos;
}

void calcSubdivPosExtra(Vertex *v, Edge *left, unsigned maxLevel)
//left vertex is extraordinary
{
	switch (left->Org()->valence)
	{
	case 3:
		calcSubdivPosExtra3(v, left, maxLevel);
		break;
	case 4:
		calcSubdivPosExtra4(v, left, maxLevel);
		break;
	default:
		calcSubdivPosExtraRest(v, left, maxLevel);
		break;
	}
}

void calcSubdivPosRegular(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;

	vec3 pos=e->Org()->subdivPos*0.5; //+1/2
	e=e->Onext();
	pos+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.125; //+1/8
	e=e->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0625; //-1/16
	e=e->Onext()->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0625; //-1/16

	e=stepAlong(left, exp2i(maxLevel - v->level + 1))->Sym(); //right
	pos+=e->Org()->subdivPos*0.5; //+1/2
	e=e->Onext();
	pos+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.125; //+1/8
	e=e->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0625; //-1/16
	e=e->Onext()->Onext();
	pos-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->subdivPos*0.0625; //-1/16
	v->subdivPos=pos;
}

void calcSubdivPos(Vertex *v, Edge *left, unsigned maxLevel)
//vertex in question and edge with origin in left neighbor (pointing towards this vertex)
{
	Edge *right=stepAlong(left, exp2i(maxLevel - v->level + 1))->Sym();
	Vertex *l=left->Org(), *r=right->Org();
	if (l->valence==6 && r->valence==6)
		calcSubdivPosRegular(v, left, maxLevel);
	else
		if (l->valence!=6 && r->valence==6)
			calcSubdivPosExtra(v, left, maxLevel);
		else if (l->valence==6 && r->valence!=6)
			calcSubdivPosExtra(v, right, maxLevel);
		else //both extraordinary
		{
			calcSubdivPosExtra(v, left, maxLevel);
			vec3 pos=v->subdivPos;
			calcSubdivPosExtra(v, right, maxLevel);
			v->subdivPos=(v->subdivPos + pos)*0.5;
		}
}

void calculateSubdivisionPositions(Cell *qe, unsigned maxLevel, unsigned subdivFromLevel)
//smooth interpolating subdivision using modified butterfly scheme
{
	for (unsigned l=subdivFromLevel+1; l<=maxLevel; l++)
	{
		//calculate subdiv positions based on previous level
		CellVertexIterator cvi(qe);
		Vertex *v;
		while ((v = cvi.next())!=0)
			if (v->level==l)
			{
				Edge *e=v->getEdge();
				Edge *e1=stepAlong(e, exp2i(maxLevel-l));
				Vertex *n=e1->Dest();
				while (n->level>=l)
				{
					e=e->Onext();
					e1=stepAlong(e, exp2i(maxLevel-l));
					n=e1->Dest();
				}
				calcSubdivPos(v, e1->Sym(), maxLevel);
			}
	}
}

void copyPos2SubdivPos(Cell *qe)
{
	CellVertexIterator cvi0(qe);
	Vertex *v0;
	while ((v0 = cvi0.next())!=0)
		v0->subdivPos=v0->pos;
}

void copySubdivPos2Pos(Cell *qe)
{
	CellVertexIterator cvi0(qe);
	Vertex *v0;
	while ((v0 = cvi0.next())!=0)
		v0->pos=v0->subdivPos;
}

void swapSubdivPos(Cell *qe)
{
	CellVertexIterator cvi0(qe);
	Vertex *v0;
	while ((v0 = cvi0.next())!=0)
		std::swap(v0->pos,v0->subdivPos);
}

void normalizeHeuristic(Cell *qe, int maxLevel, int howManyFinestLevels)
{
	Vertex *v;
	//objWriteCell(qe, "D:\\Temp\\hierarchyP.obj");

	for (int l=maxLevel-1; l>maxLevel-howManyFinestLevels; l--)
	{
		double weight=0.5/pow(3.0,int(maxLevel-l));
		CellVertexIterator cvi(qe);
		while ((v = cvi.next())!=0)
			if (v->level==l)
			{
				Edge *e=v->getEdge();
				vec3 update(0,0,0);
				do
				{
					Vertex *n=stepAlong(e, exp2i(maxLevel-l-1))->Dest();
                    update+=(n->pos - n->subdivPos);
					e=e->Onext();
				} while (e!=v->getEdge());
                v->pos=v->subdivPos+(v->pos - v->subdivPos)*weight+(update/v->valence)*(1-weight);
			}
	}
	
    //all coarse levels in one go
    double weight=0.5/3;
	CellVertexIterator cvi(qe);
	while ((v = cvi.next())!=0)
		if (v->level<=maxLevel-howManyFinestLevels)
		{
			Edge *e=v->getEdge();
			vec3 update(0,0,0);
			do
			{
				Vertex *n=e->Dest();
                update+=(n->pos - n->subdivPos);
				e=e->Onext();
			} while (e!=v->getEdge());
            v->pos=v->subdivPos+(v->pos - v->subdivPos)*weight+(update/v->valence)*(1-weight);
		}

	//objWriteCell(qe, "D:\\Temp\\hierarchyN.obj");
	copyPos2SubdivPos(qe);
	calculateSubdivisionPositions(qe, maxLevel, maxLevel-howManyFinestLevels);
	copySubdivPos2Pos(qe);
	//objWriteCell(qe, "D:\\Temp\\hierarchyS.obj");
}

void calcWeightExtra3(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	v->weights=e->Org()->weights*0.75; //+9/12
	v->weights+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.4166666666666666666; //+5/12
	e=e->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0833333333333333333; //-1/12
	e=e->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0833333333333333333; //-1/12
}

void calcWeightExtra4(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	v->weights=e->Org()->weights*0.75; //+6/8
	v->weights+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.375; //+3/8
	e=e->Onext()->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.125; //-1/8
}

void calcWeightExtraRest(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;
	v->weights=e->Org()->weights*0.75;
	double tc;
	unsigned valence=left->Org()->valence;
	for (unsigned i=0; i<valence; i++)
	{
		tc=1.0/valence*(0.25+cos(2*i*M_PI/valence)+0.5*cos(4*i*M_PI/valence));
		v->weights+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*tc;
		e=e->Onext();
	}
}

void calcWeightExtra(Vertex *v, Edge *left, unsigned maxLevel)
//left vertex is extraordinary
{
	switch (left->Org()->valence)
	{
	case 3:
		calcWeightExtra3(v, left, maxLevel);
		break;
	case 4:
		calcWeightExtra4(v, left, maxLevel);
		break;
	default:
		calcWeightExtraRest(v, left, maxLevel);
		break;
	}
}

void calcWeightRegular(Vertex *v, Edge *left, unsigned maxLevel)
{
	Edge *e=left;

	v->weights=e->Org()->weights*0.5; //+1/2
	e=e->Onext();
	v->weights+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.125; //+1/8
	e=e->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0625; //-1/16
	e=e->Onext()->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0625; //-1/16

	e=stepAlong(left, exp2i(maxLevel - v->level + 1))->Sym(); //right
	v->weights+=e->Org()->weights*0.5; //+1/2
	e=e->Onext();
	v->weights+=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.125; //+1/8
	e=e->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0625; //-1/16
	e=e->Onext()->Onext();
	v->weights-=stepAlong(e, exp2i(maxLevel - v->level + 1))->Dest()->weights*0.0625; //-1/16
}

void calcWeight(Vertex *v, Edge *left, unsigned maxLevel)
//vertex in question and edge with origin in left neighbor (pointing towards this vertex)
{
	Edge *right=stepAlong(left, exp2i(maxLevel - v->level + 1))->Sym();
	Vertex *l=left->Org(), *r=right->Org();
	if (l->valence==6 && r->valence==6)
		calcWeightRegular(v, left, maxLevel);
	else
		if (l->valence!=6 && r->valence==6)
			calcWeightExtra(v, left, maxLevel);
		else if (l->valence==6 && r->valence!=6)
			calcWeightExtra(v, right, maxLevel);
		else //both extraordinary
		{
			calcWeightExtra(v, left, maxLevel);
			Eigen::VectorXd temp=v->weights;
			calcWeightExtra(v, right, maxLevel);
			v->weights=(v->weights + temp)*0.5;
		}
}

Eigen::MatrixXd createWeightMatrix(Cell *qe, unsigned maxLevel, unsigned subdivFromLevel, unsigned &freeCount)
{
	Vertex *v;
	freeCount=0;
	{
		CellVertexIterator cvi0(qe);
		while ((v = cvi0.next())!=0)
			if (v->level<=subdivFromLevel)
				freeCount++;
	} //we now have count of free vertices

	CellVertexIterator cvi0(qe);
	while ((v = cvi0.next())!=0)
		if (v->level<=subdivFromLevel)
		{
			v->weights.resize(freeCount);
			v->weights.setZero();
			v->weights[v->realIndex()]=1;
		}

	for (unsigned l=subdivFromLevel+1; l<=maxLevel; l++)
	{
		CellVertexIterator cvi0(qe);
		while ((v = cvi0.next())!=0)
			if (v->level==l)
			{
				v->weights.resize(freeCount);
				v->weights.setZero();
				Edge *e=v->getEdge();
				Edge *e1=stepAlong(e, exp2i(maxLevel-l));
				Vertex *n=e1->Dest();
				while (n->level>=l)
				{
					e=e->Onext();
					e1=stepAlong(e, exp2i(maxLevel-l));
					n=e1->Dest();
				}
				calcWeight(v, e1->Sym(), maxLevel);
			}
	}
	
	Eigen::MatrixXd W(freeCount, qe->countVertices()); //freeCount rows, totalCount columns
	CellVertexIterator cvi(qe);
	while ((v = cvi.next())!=0)
	{
		W.col(v->realIndex())=v->weights.transpose(); //vectors are columns, not rows (mathematical standard)
	}
	return W;
}