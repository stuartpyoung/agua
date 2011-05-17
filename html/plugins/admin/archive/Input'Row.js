dojo.provide( "plugins.admin.ParameterRow");


dojo.declare( "plugins.admin.ParameterRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/parameterrow.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Apps WIDGET
	parentWidget : null,

	// FORM INPUTS AND TYPES (word|phrase)
	formInputs : {
        name: "word",
        argument: "word",
        type: "word",
        category: "word",
        value: "word",
        discretion: "word",
        description: "phrase",
        format: "word",
        args: "word",
        params: "phrase",
        paramFunction: "phrase"
	},


	constructor : function(args)
	{

		this.parentWidget = args.parentWidget;

		//this.inherited(arguments);
	},


	postCreate : function()
	{

		this.startup();
	},

	startup : function ()
	{

		this.inherited(arguments);

		// CONNECT TOGGLE EVENT
		var parameterRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {
			parameterRowObject.toggle();
		});

		// ADD 'ONCLICK' EDIT VALUE LISTENERS
		var parameterRowObject = this;
		var array = [ "argument", "category", "value", "description", "format", "args", "params", "paramFunction" ];
		for ( var i in array )
		{
			dojo.connect(this[array[i]], "onclick", function(event)
				{
					parameterRowObject.parentWidget.editParameterRow(parameterRowObject, event.target);
					event.stopPropagation(); //Stop Event Bubbling
				}
			);
		}

		// ADD 'ONCHANGE' COMBO BOX LISTENERS
		var parameterRowObject = this;
		var array = [ "type", "discretion", "paramtype" ];
		for ( var i in array )
		{
			dojo.connect(this[array[i]], "onchange", function(event)
				{
					var inputs = parameterRowObject.parentWidget.getEditedInputs(parameterRowObject);
					parameterRowObject.parentWidget.saveParameter(inputs, null);
					event.stopPropagation(); //Stop Event Bubbling
				}
			);
		}
	},

	// TOGGLE HIDDEN NODES
	toggle : function ()
	{

		var array = [ "category", "argument", "type", "value", "description", "format", "args", "params", "paramFunction" ];
		for ( var i in array )
		{
			if ( this[array[i]].style.display == 'table-cell' ) this[array[i]].style.display='none';
			else this[array[i]].style.display = 'table-cell';
		}
	}

});

