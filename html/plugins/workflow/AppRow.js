
dojo.provide( "plugins.workflow.AppRow");


dojo.declare( "plugins.workflow.AppRow",
	[ dijit._Widget, dijit._Templated ],
{
	/////}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/approw.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT plugins.workflow.Apps WIDGET
parentWidget : null,

// CORE WORKFLOW OBJECTS
core : null,


// APPLICATION OBJECT
application : null,

constructor : function(args)
{

	this.core = args.core;
	this.parentWidget = args.parentWidget;

	this.application = new Object;
	for ( var key in args )
	{
		if ( key != "parentWidget" )
		{
			this.application[key] = args[key];
		}
	}

	//this.inherited(arguments);
},

// RETURN A COPY OF this.application
getApplication : function ()
{
	return dojo.clone(this.application);
},

// SET this.application TO THE SUPPLIED APPLICATION OBJECT
setApplication : function (application)
{
	this.application = application;

	return this.application;
},



postCreate : function()
{

	this.startup();
},

startup : function ()
{

	this.inherited(arguments);


	// HACK:
	//
	// SET parentWidget TO this.name FOR RETRIEVAL OF this.application
	// WHEN MENU IS CLICKED
	//
	// REM: remove ONCLICK BUBBLES ON appRow.name NODE RATHER THAN ON node. 
	// I.E., CONTRARY TO DESIRED, this.name IS THE TARGET INSTEAD OF THE node.
	//
	// ALSO ADDED node.parentWidget = appRow IN Stages.setDropTarget()

	this.name.parentWidget = this;

	// CONNECT TOGGLE EVENT
	var appRowObject = this;
	dojo.connect( this.name, "onclick", function(event) {
		event.stopPropagation();
		appRowObject.toggle();
	});

},

toggle : function ()
{

	var array = [ "description" ];
	for ( var i in array )
	{
		if( this[array[i]].style )
		{
			if ( this[array[i]].style.display == 'table-cell' ) this[array[i]].style.display='none';
			else this[array[i]].style.display = 'table-cell';
		}
	}
}




});


