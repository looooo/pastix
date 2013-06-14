/************************************************************/
/**                                                        **/
/**   NAME       : ps.h                                    **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                write postscript with draws             **/
/**                of symbol matrix and elimination tree   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/

#ifndef PS_H
#define static
#endif

FILE    *       ps_open(char *);
void            ps_close(FILE *);
void            ps_write_matrix(SymbolMatrix *, FILE *, PASTIX_INT *);
void            ps_write_tree(const CostMatrix *, const EliminTree *, FILE *, PASTIX_INT *);
static double    ps_rec_write_tree(PASTIX_INT , const CostMatrix *, const EliminTree *, FILE *,
			    void (*ps_draw_node)(FILE *,PASTIX_INT ,const CostMatrix *, const EliminTree *,double ,double ,double )
				  );
static void     ps_draw_node_num(FILE *, PASTIX_INT , const CostMatrix *,const EliminTree *, double , double , double );
void            ps_write_tree_owner(PASTIX_INT *,const CostMatrix *, const EliminTree *, FILE *, PASTIX_INT *);
static double    ps_rec_write_tree_owner(PASTIX_INT , PASTIX_INT *, const CostMatrix *, const EliminTree *, FILE *,
			    void (*ps_draw_node)(FILE *,PASTIX_INT, PASTIX_INT, const CostMatrix *, const EliminTree *,double ,double ,double )
				  );
static void     ps_draw_node_owner(FILE *, PASTIX_INT, PASTIX_INT, const CostMatrix *,const EliminTree *, double , double , double );




#undef static
