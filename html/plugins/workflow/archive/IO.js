dojo.provide( "plugins.workflow.IO");

/* CLASS SUMMARY: LOGICALLY CHAIN TOGETHER INPUTS/OUTPUTS OF WORKFLOW STAGES

  detailed methodology for Workflow IO

	METHOD:		workflowIO

	FOR EACH APPLICATION IN stages

 		1. SET INPUTS DEPENDENT ON THE OUTPUTS OF THE PREVIOUS APPLICATION

 		2. SET RESOURCES DEPENDENT ON INPUTS

 		3. SET OUTPUTS DEPENDENT ON INPUTS AND RESOURCES

 (NB: REPEAT STEP 3 IF USER CHANGES ARGUMENTS MANUALLY)



 THE REASON FOR EXTRACTING THE LIST OF OUTPUT FILES FROM THE
 APPLICATION LIKE THIS IS THAT THE OUTPUT FILES MAY NOT BE
 MENTIONED IN THE ARGUMENTS SO WE HAVE TO INFER THE OUTPUT FILE
 NAMES FROM THE APPLICATIONS ARGUMENT LIST USING SOME ADDITIONAL LOGIC
 CONTAINED IN THE 'outputs' ENTRY FOR EACH APPLICATION. E.G.:
  "outputs":
	{
		{
			'name': "acefile",
			'format': 'ace', 'type': 'file',
			'arguments': [ 'arguments.inputfile.value' ],
			'javascript': "function(inputfile) { var acefile = inputfile; acefile = acefile.replace(/\\.txt/, \".ace\"); return acefile; }"
		}
	}


	"resources":
	{
		"sortedfile":
		{
			'format': 'sorted',
			'type': 'file',
			'args': [ 'outputs.sortedfile.value' ],
			'javascript': "new Function('inputfile', 'var acefile = inputfile; acefile = acefile.replace(/\\.txt/, \".ace\"); return acefile;')"
		}
	},


	'arguments':
	{
		'type':
		{
			'option': '--type',
			'value': 'single',
			'type': 'parameter',
			'description': "Either 'single' or 'paired' reads",
			'required' : true
		},

		'rundir':
		{
			'option': '--rundir',
			'value': 'Project1/Workflow3-indels',
			'type': 'file',
			'format': 'eland-sorted',
			'description': '',
			'discretion' : 'required'
		}
	}
*/


// INHERITS
dojo.require("plugins.core.Common");

dojo.declare( "plugins.workflow.IO",
	[ plugins.core.Common ],
{
	// PARENT WIDGET
	parentWidget : null,

	////}}

constructor : function(args) {
	// GET INFO FROM ARGS
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

	// SET DRAG APP - LIST OF APPS
	this.setMenu();
},

chainInput : function () {
// GET THE INPUT FOR AN APPLICATION DEPENDING ON THE OUTPUTS
// AND RESOURCES OF PRECEDING APPLICATIONS IN THE WORKFLOW

},

reload : function (stages) {
// GET THE OUTPUT 

	for ( var stageIndex = 0; stageIndex < 1; stageIndex++ )
	{
		var stage = stages[stageIndex];
		var previousApplication = stages[stageIndex - 1];

		// 1. SET INPUTS DEPENDENT ON THE OUTPUTS OF THE PREVIOUS APPLICATION
		if ( stageIndex > 0 )
		{
			this.getValues(stage.inputs, previousApplication);
		} // if ( stageIndex > 0 )

		// 2. SET ARGUMENTS DEPENDENT ON INPUTS
		this.getValues(stage.arguments, stage);

		// 3. SET OUTPUTS DEPENDENT ON ARGUMENTS
		this.getValues(stage.outputs, stage);

	} // foreach stages
},

getValues : function(resources, stage) {

	// GET resources
	for ( var resource in resources )
	{

		var functionArguments = resources[resource].args;

		var parameterValues = this.parameterValues(stage, functionArguments);


		//if ( resources[resource].function )
		//{
		//}

		if ( parameterValues )
		{

//
//				// CLEAN resources[resource].paramFunction
//				resources[resource].paramFunction = this.jsonSafe(resources[resource].paramFunction, 'fromJson');
//

			// CLEAN resources[resource].paramFunction
			//resources[resource].paramFunction = resources[resource].paramFunction.replace(/&quot;/g, "'");


			var value = this.value(parameterValues, resources[resource].params, resources[resource].paramFunction);



			if ( value && value != '' && value != "undefined" && value[0] != "undefined" )
			{
				resources[resource]["value"] = value;
			}
			else
			{
			}
		}
		else
		{
		}




	}
},

parameterValues : function (stage, functionArguments) {
/* METHOD: parameterValues

 EXTRACT ARGUMENT VALUES FROM APPLICATION OBJECT FOR EACH ARGUMENT
 USED AS FOLLOWS:

		CURRENT stage AND CURRENT APPLICATION functionArguments 
			-> GENERATE OUTPUTS FOR CURRENT APPLICATION

		PREVIOUS stage AND CURRENT APPLICATION functionArguments 
			-> GENERATE INPUTS FOR CURRENT APPLICATION
*/


	if ( ! functionArguments )
	{
		return;
	}

	// CONVERT functionArguments TO ARRAY


	// CONVERT ARRAY/HASH ARGUMENTS FROM TEXT INTO OBJECTS
	if ( typeof functionArguments != "ARRAY" && typeof functionArguments != "object")
	{
		try {
			functionArguments = eval(functionArguments);	
		}
		catch (e)
		{
			functionArguments = '';
		}

	}

	// RETURN AN ARRAY
	var parameterValues = new Array;
	for ( var index = 0; index < functionArguments.length; index++ )
	{

		// EXTRACT ARGUMENT VALUE FROM APPLICATION OBJECT
		var array = functionArguments[index].split("\.");
		var value = stage;
		if ( value == null || value == ''  || value == "undefined" )
		{
			parameterValues.push('');
		}
		else
		{
			for ( var itemIndex = 0; itemIndex < array.length; itemIndex++ )
			{
				if ( value == "undefined" || value == null || value[array[itemIndex]] == null )
				{
					parameterValues.push('');
					break;
				}
				else
				{
					var tempValue = value[array[itemIndex]];
					value = tempValue;
				}
			}	
			parameterValues.push(value);
		}
	}


	return parameterValues;
},

value : function (parameterValues, params, newFunction) {

	newFunction = this.jsonSafe(newFunction, 'fromJson');

	var value;
	if ( ! params )
	{
		if ( parameterValues )
		{
			value = parameterValues[0];
		}
		else
		{
			value = '';
		}
	}
	else
	{

		var argumentFunction = eval("new Function(\"" + params + "\", \"" + newFunction + "\")");



		//var argumentString  = "'s_6_150' " + ", 'Project1/Workflow3-indels'";
//        var argumentString = '';
//        for ( var index in parameterValues )
//        {
//            argumentString += "'" + parameterValues[index] + "',";
//        }
//		argumentString = argumentString.replace(/,$/, '');
		//value = argumentFunction(argumentString);

		//var array = new Array;
		//array.push(parameterValues[0]);
		//array.push(parameterValues[1]);
		//value = argumentFunction(array);


		// HACK TO DEAL WITH UP TO FIVE ARGUMENTS
		value = argumentFunction(parameterValues[0], parameterValues[1], parameterValues[2], parameterValues[3], parameterValues[4]);

		// DIDN'T WORK
		//value = argumentFunction(parameterValues.split(","));
		//value = argumentFunction(parameterValues[0], parameterValues[1]);

		// WORKED
		//value = argumentFunction('s_6_150', 'Project1/Workflow3-indels');
	}

	return value;
}

});
