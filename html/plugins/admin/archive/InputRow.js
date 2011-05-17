dojo.provide( "plugins.admin.InputRow");


dojo.declare( "plugins.admin.InputRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/inputrow.html"),

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
		var inputRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {
			inputRowObject.toggle();
		});

		// ADD 'ONCLICK' EDIT VALUE LISTENERS
		var inputRowObject = this;
		var array = [ "argument", "category", "value", "description", "format", "args", "params", "paramFunction" ];
		for ( var i in array )
		{
			dojo.connect(this[array[i]], "onclick", function(event)
				{
					inputRowObject.parentWidget.editInputRow(inputRowObject, event.target);
					event.stopPropagation(); //Stop Event Bubbling
				}
			);
		}

		// ADD 'ONCHANGE' COMBO BOX LISTENERS
		var inputRowObject = this;
		var array = [ "type", "discretion", "paramtype" ];
		for ( var i in array )
		{
			dojo.connect(this[array[i]], "onchange", function(event)
				{
					var inputs = inputRowObject.parentWidget.getEditedInputs(inputRowObject);
					inputRowObject.parentWidget.saveInput(inputs, null);
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

