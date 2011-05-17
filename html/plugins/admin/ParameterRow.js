dojo.provide( "plugins.admin.ParameterRow");


dojo.declare( "plugins.admin.ParameterRow",
	[ dijit._Widget, dijit._Templated ],
{
	////}}
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/parameterrow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT plugins.admin.Apps WIDGET
parentWidget : null,

constructor : function(args) {

	this.parentWidget = args.parentWidget;
	this.lockedValue = args.locked;
},

postCreate : function(args) {
	//
	//this.parameter = parameter;
	this.startup();
},

startup : function () {

	this.inherited(arguments);

	// CONNECT TOGGLE EVENT
	var parameterRowObject = this;
	dojo.connect( this.name, "onclick", function(event) {
		parameterRowObject.toggle();
	});

	// SET LOCKED CLASS
	if ( this.lockedValue == 1 )	{
		dojo.removeClass(this.locked,'unlocked');
		dojo.addClass(this.locked,'locked');
	}
	else	{
		dojo.removeClass(this.locked,'locked');
		dojo.addClass(this.locked,'unlocked');
	}

	// ADD 'ONCLICK' EDIT VALUE LISTENERS
	var parameterRowObject = this;
	var onclickArray = [ "argument", "category", "value", "description", "format", "args", "params", "paramFunction" ];
	for ( var i in onclickArray )
	{
		dojo.connect(this[onclickArray[i]], "onclick", function(event)
			{
				parameterRowObject.parentWidget.editParameterRow(parameterRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling
			}
		);
	}

	// ADD 'ONCHANGE' COMBO BOX LISTENERS
	var parameterRowObject = this;
	var onchangeArray = [ "valuetype", "ordinal", "discretion", "paramtype" ];
	for ( var i in onchangeArray )
	{
		dojo.connect(this[onchangeArray[i]], "onchange", function(event)
			{
				var inputs = parameterRowObject.parentWidget.getEditedInputs(parameterRowObject);
				parameterRowObject.parentWidget.saveParameter(inputs, null);
				event.stopPropagation(); //Stop Event Bubbling
			}
		);
	}
},

toggle : function () {
// TOGGLE HIDDEN NODES

	var array = [ "category", "argument", "valuetype", "value", "ordinal", "description", "format", "args", "params", "paramFunction" ];
	for ( var i in array )
	{
		if ( this[array[i]].style.display == 'table-cell' ) this[array[i]].style.display='none';
		else this[array[i]].style.display = 'table-cell';
	}
},

toggleLock : function (event) {

	if ( dojo.hasClass(this.locked, 'locked')	){
		dojo.removeClass(this.locked, 'locked');
		dojo.addClass(this.locked, 'unlocked');
		Agua.toastMessage("ParameterRow has been unlocked. Users can change this parameter.")
	}	
	else {
		dojo.removeClass(this.locked, 'unlocked');
		dojo.addClass(this.locked, 'locked');

		Agua.toastMessage("ParameterRow has been locked. Users will not be able to change this paramter.")
	}	

	var inputs = this.parentWidget.getEditedInputs(this);
	this.parentWidget.saveParameter(inputs, null);
	event.stopPropagation(); //Stop Event Bubbling
},

formInputs : {
// FORM INPUTS AND TYPES (word|phrase)
	locked : "",
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
}




});
