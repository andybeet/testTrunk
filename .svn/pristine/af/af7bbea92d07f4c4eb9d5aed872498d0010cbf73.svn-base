/**
 *	\file The C header file containing the atPhysicsStructure definition.
 *
 *
 *
 *
 */
typedef struct {
	/**
	 * @name The sourceSink Files. These are opened in sourceSink
	 * and left open. They are closed at the end of the model run in the function
	 * freePhysicsStruct().
	 */
	//@{
	FILE *inpfp;
	double *totinp;
	//@}

	/**
	 *
	 */
	double **masstosed;

	/* Transport files*/
	FILE *expfp;

}atPhysicsStructure;
