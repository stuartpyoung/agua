dojo.provide( "plugins.workflow.IO");

// SET THE DEFAULT CHAINED VALUES FOR INPUTS AND OUTPUTS FOR AN
// APPLICATION BASED ON THOSE OF THE PREVIOUS APPLICATIONS

// INHERITS
dojo.require("plugins.core.Common");

dojo.declare( "plugins.workflow.IO",
	[ plugins.core.Common ],
{
	////}

// PARENT WIDGET
parentWidget : null,

// CORE WORKFLOW OBJECTS
core : null,

constructor : function(args) {
	// GET INFO FROM ARGS
	this.core = args.core;
	this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 
},

// GET INPUTS FOR AN APPLICATION, USUALLY BASED ON:
//      THE OUTPUTS AND RESOURCES OF THE PREVIOUS APPLICATION
chainInputs: function (application, force) {

	// GET THE input STAGE PARAMETERS FOR THIS STAGE        
	var stageParameters = Agua.getStageParameters(application);
	var inputParameters = this.filterByKeyValues(stageParameters, ["paramtype"], ["input"]);

	// GET THE STAGE PARAMETERS FOR THE PRECEDING STAGE
	var appnumber = application.appnumber;
	if ( appnumber == null )
	{
		return;
	}

	var stages = Agua.getStagesByWorkflow(application.project,application.workflow);

	if ( stages == null )
	{
		return;
	}

	var previousStage = stages[appnumber - 2];
	if ( previousStage == null )
	{
		return;
	}

	var previousStageParameters = Agua.getStageParameters(previousStage);
	if ( previousStageParameters == null )
	{
		return;
	}

	for ( var i = 0; i < inputParameters.length; i++ )
	{
		this.chainStageParameter(inputParameters[i], previousStageParameters, force);
	}

},  //  chainInputs


// GET RESOURCES FOR AN APPLICATION, USUALLY BASED ON:
//      ITS INPUTS
chainResources: function (application, force) {


},


// GET OUTPUTS FOR AN APPLICATION, USUALLY BASED ON:
//      ITS OWN INPUTS AND RESOURCES
chainOutputs : function (application, force) {

	// GET THE input STAGE PARAMETERS FOR THIS STAGE        
	var stageParameters = Agua.getStageParameters(application);

	// GET OUTPUT STAGE PARAMETERS ONLY
	var outputParameters = dojo.clone(stageParameters);
	outputParameters = this.filterByKeyValues(outputParameters, ["paramtype"], ["output"]);

	for ( var i = 0; i < outputParameters.length; i++ )
	{
		this.chainStageParameter(outputParameters[i], stageParameters, force);
	}
},

chainStageParameter : function (stageParameter, sourceParameters, force) {

	// RETURN IF args IS NULL
	if ( stageParameter.args == null || stageParameter.args == '' )
	{
		return;
	}

	var valuesArray = this.getValuesArray(stageParameter, dojo.clone(sourceParameters));
	var value = this.getValue(valuesArray, stageParameter.params, stageParameter.paramFunction);

	// REPLACE EMPTY STAGE PARAMETER VALUE IN Agua AND REMOTE DATABASE 
	// OR OVERWRITE EXISTING VALUE IF force IS TRUE
	if ( (force == true && value != null)
		|| (value != null && stageParameter.value == '') )
	{
		stageParameter.value = value;
		stageParameter.chained = 1;
		Agua._removeStageParameter(stageParameter);
		Agua.addStageParameter(stageParameter);
	}
},


// 		1. SET INPUTS DEPENDENT ON THE OUTPUTS OF THE PREVIOUS STAGE
// 		2. SET RESOURCES DEPENDENT ON INPUTS OF THIS STAGE
// 		3. SET OUTPUTS DEPENDENT ON INPUTS AND RESOURCES OF THIS STAGE
//      (NB: REPEAT STEP 3 IF USER CHANGES RESOURCES MANUALLY)
chainStage : function (application, force) {

	// SET INPUTS
	// RETURN IF NO PRECEDING APPLICATIONS
	if ( application.number != 1 )
	{
		this.chainInputs(application, force);
	}

	//// SET RESOURCES
	//application = this.chainResources(application.arguments, application);
	//

	// SET OUTPUTS 
	this.chainOutputs(application, force);


	return;
},



// GET THE ARRAY OF VALUES FOR THE args TO BE INPUT INTO THE paramFunction
getValuesArray : function(parameter, sourceParameters) {

	// SANITY CHECK		
	if ( parameter.args == null || sourceParameters == null )
	{
		return;
	}

	// CONVERT args TO argsArray
	var argsArray = parameter.args.split(/,/);

	// RETURN AN ARRAY
	var valuesArray = new Array;
	for ( var i = 0; i < argsArray.length; i++ )
	{

		// EXTRACT ARGUMENT VALUE FROM APPLICATION OBJECT
		var array = argsArray[i].split("\.");

		var paramtypeToken = array[0];
		var nameToken = array[1];
		var valueToken = array[2];

		// FILTER BY PARAMTYPE TYPE (input|resource|output)
		sourceParameters = this.filterByKeyValues(sourceParameters, ["paramtype"], [paramtypeToken]);

		// IF SOURCE PARAMETERS NOT PRESENT, SET VALUE TO ''
		if ( sourceParameters == null || sourceParameters.length == 0 )
		{
			valuesArray.push('');
			continue;
		}

		// IF THE ARRAY OF SOURCE PARAMETERS WAS WANTED, STORE IT AND NEXT
		if ( nameToken == null || nameToken == '' )
		{
			valuesArray.push(sourceParameters);
			continue;
		}

		// GET THE SOURCE PARAMETER
		var sourceParameter = this.filterByKeyValues(sourceParameters, ["name"], [nameToken])[0];

		// IF SOURCE PARAMETER NOT DEFINED, SET VALUE TO ''
		if ( sourceParameter == null || sourceParameter == '' )
		{
			valuesArray.push('');
			continue;
		}

		// IF THE WHOLE PARAMETER HASH WAS WANTED, STORE IT AND NEXT
		if ( valueToken == null )
		{
			valuesArray.push(sourceParameter);
			continue;
		}

		// OTHERWISE, GET THE PARAMETER'S VALUE
		var value = sourceParameter.value;
		valuesArray.push(value);
	}


	return valuesArray;
},


// GET THE PARAMETER VALUE USING THE paramFunction IF PRESENT.
// OTHERWISE, RETURN A SCALAR OF THE FIRST ENTRY IN parameterValues
getValue : function (parameterValues, params, paramFunction) {

	var value;
	if ( params == null || params == ''
		|| paramFunction == null || paramFunction == '' )
	{
		return parameterValues[0];
	}

	// SET THE FUNCTION
	var paramsArray = params.slice(',');
	var argumentFunction = new Function( paramsArray , paramFunction );

	// RUN THE FUNCTION WITH THE INPUT PARAMETER VALUES
	// HACK TO INPUT PARAMETER VALUES
	value = argumentFunction(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3], parameterValues[4]);
	// THIS DOESN'T WORK
	//value = argumentFunction(parameterValues);

	return value;
}

});

/*

        stageParameters = [ 	
	{
        "owner": "admin",
		 "value": "",
		 "args": "",
		 "name": "inputfile",
		 "description": "",
		 "discretion": "required",
		 "params": "",
		 "appname": "image2eland.pl",
		 "apptype": "pipeline",
		 "argument": "--inputfile",
		 "format": "",
		 "paramtype": "input",
		 "category": "inputfile",
		 "valuetype": "file",
		 "paramFunction": "",
		 "project": "aaa",
		 "workflow": "test",
		 "appnumber": "2"
    },

	{
        "owner": "admin",
		 "value": "",
		 "args": "",
		 "name": "outputdir",
		 "description": "",
		 "discretion": "required",
		 "params": "",
		 "appname": "image2eland.pl",
		 "apptype": "pipeline",
		 "argument": "--outputdir",
		 "format": "",
		 "paramtype": "input",
		 "category": "outputdir",
		 "valuetype": "directory",
		 "paramFunction": "",
		 "project": "aaa",
		 "workflow": "test",
		 "appnumber": "2"
    },

	{
        "owner": "admin",
		 "value": "",
		 "args": "input.acefile.value",
		 "params": "",
		 "paramFunction": "var output",

		 "name": "outputfile",
		 "description": "",
		 "discretion": "required",
		 "appname": "image2eland.pl",
		 "apptype": "pipeline",
		 "argument": "--outputfile",
		 "format": "",
		 "paramtype": "output",
		 "category": "outputfile",
		 "valuetype": "file",
		 "project": "aaa",
		 "workflow": "test",
		 "appnumber": "2"
    },

	{
        "owner": "admin",
		 "value": "",
		 "args": "input.acefile.value",
		 "params": "",
		 "paramFunction": "var output",

		 "name": "outputfile",
		 "description": "",
		 "discretion": "required",
		 "appname": "image2eland.pl",
		 "apptype": "pipeline",
		 "argument": "--outputfile",
		 "format": "",
		 "paramtype": "output",
		 "category": "outputfile",
		 "valuetype": "file",
		 "project": "aaa",
		 "workflow": "test",
		 "appnumber": "2"
    }
];


// E.G., THE NAME OF THE OUTPUT FILE FOR AN APPLICATION OFTEN
// DEPENDS ON THE INPUTS FOR THE APPLICATION. IN THIS CASE, THE
// args, params AND paramFunction ENTRIES PROVIDE THE MEANS
// TO INFER THE NAME OF THE OUTPUT FILE:

	{
		{
			name: "acefile",
          paramtype: "output",
			format: "ace",
          type: "file",
			args: [ 'arguments.inputfile.value' ],
          params : [inputfile],
			paramFunction: "var acefile = inputfile; acefile = acefile.replace(/\\.txt/, \".ace\"); return acefile;"
		},

		{
			name: "sortedfile":,
          format: 'sorted',
          paramtype: "output",
			type: "file",
			args: [ 'outputs.sortedfile.value' ],
			paramFunction: "new Function('inputfile', 'var acefile = inputfile; acefile = acefile.replace(/\\.txt/, \".ace\"); return acefile;')"
		},

		{
			name: 'type:
		    option: '--type',
			value: 'single',
			type: 'parameter',
			description: "Either 'single' or 'paired' reads",
			'required: true
		},

		{
			option: '--rundir',
      	name: 'rundir',
			value: 'Project1/Workflow3-indels',
			type: 'file',
			format: 'eland-sorted',
			description: '',
			discretion: 'required'
		}
	}

*/
