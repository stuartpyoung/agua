dojo.provide( "plugins.admin.ClusterRow");


dojo.declare( "plugins.admin.ClusterRow",
	[ dijit._Widget, dijit._Templated ],
{
	////}}
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/clusterrow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT plugins.admin.Apps WIDGET
parentWidget : null,

// FORM INPUTS AND TYPES (word|phrase)
formInputs : {
	clusterName	:	"word",
	minNodes	:	"word",
	maxNodes	:	"word",
	typeCombo	:	"word",
	amiId		:	"word",
	description	:	"phrase",
	notes		:	"phrase"
},


constructor : function(args) {
	this.parentWidget = args.parentWidget;
	this.lockedValue = args.locked;
},


postCreate : function() {

	this.startup();
},

startup : function () {

	this.inherited(arguments);

	//// ADD 'ONCLICK' EDIT VALUE LISTENERS
	//var thisObject = this;
	//var onclickArray = [ "minNodes", "maxNodes", "amiId", "description", "notes" ];
	//for ( var i in onclickArray )
	//{
	//	dojo.connect(this[onclickArray[i]], "onclick", function(event)
	//		{
	//			thisObject.parentWidget.editClusterRow(thisObject, event.target);
	//			event.stopPropagation(); //Stop Event Bubbling
	//		}
	//	);
	//}

	// ADD 'ONCHANGE' COMBO BOX LISTENERS
	var thisObject = this;
	var onchangeArray = [ "typeCombo" ];
	for ( var i in onchangeArray )
	{
		dojo.connect(this[onchangeArray[i]], "onchange", function(event)
			{
				var inputs = thisObject.parentWidget.getEditedInputs(thisObject);
				thisObject.parentWidget.saveCluster(inputs, null);
				event.stopPropagation(); //Stop Event Bubbling
			}
		);
	}
},

saveCluster : function () {

	this.checkNodeNumbers();

	var inputs = this.parentWidget.getEditedInputs(this);

	this.parentWidget.saveCluster(inputs, null);
},

checkNodeNumbers : function() {

	if (this.minNodes.get('value') > this.maxNodes.get('value') )
	{
		this.minNodes.set('value', this.maxNodes.get('value'));
	}
},



editCluster : function (event) {

	this.parentWidget.editClusterRow(this, event.target);
	event.stopPropagation(); //Stop Event Bubbling
},

checkEnterNodes : function (event) {

	if (event.keyCode == dojo.keys.ENTER)
	{
		document.body.focus();

		this.checkNodeNumbers();

		var inputs = this.parentWidget.getEditedInputs(this);
		this.saveCluster(inputs, null);

		dojo.stopEvent(event);
	}
},



// TOGGLE HIDDEN NODES
toggle : function () {

	var array = [ "typeComboContainer", "amiId", "description", "notes" ];
	for ( var i in array )
	{
		if ( this[array[i]].style.display == 'table-cell' )
			this[array[i]].style.display='none';
		else
			this[array[i]].style.display = 'table-cell';
	}
}



});

