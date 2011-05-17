dojo.provide( "plugins.core.SelectiveDialog");

/* CLASS SUMMARY: AN INTERACTIVE DIALOG WITH AN OPTIONAL COMBOBOX 

	AND OPTIONAL CHECKBOX.

	LIKE ITS INHERITED CLASS, InteractiveDialog, SelectiveDialog

	WAITS UNTIL THE enterCallback METHOD CLOSES IT, SO THE 

	enterCallback METHOD CAN VALIDATE THE INPUT AND CLOSE

	THE DIALOG WHEN THE CORRECT INPUT IS PRESENT.

*/

// HAS A
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");
dojo.require("plugins.core.SelectiveDialog");


dojo.declare( "plugins.core.SelectiveDialog",
	[ plugins.core.InteractiveDialog ],
{
	//////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "core/templates/selectivedialog.html"),

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins.core") + "/css/selectivedialog.css" ],

constructor : function(args) {

	this.title 				=	args.title;
	this.message 			=	args.message;
	this.parentWidget 		=	args.parentWidget;
	this.enterCallback 		=	args.enterCallback;
	this.cancelCallback 	=	args.cancelCallback;

	// LOAD CSS
	this.loadCSS();
},

// ADD CSS NAMESPACE CLASS
setNamespaceClass : function () {
	dojo.addClass(this.dialog.containerNode, "selectiveDialog");
	dojo.addClass(this.dialog.titleNode, "selectiveDialog");
	dojo.addClass(this.dialog.closeButtonNode, "selectiveDialog");	
},

postCreate : function() {
	this.startup();
},

startup : function () {

	this.inherited(arguments);

	// SET UP DIALOG
	this.setDialogue();

	// ADD CSS NAMESPACE CLASS
	dojo.addClass(this.dialog.containerNode, "inputDialog");
	dojo.addClass(this.dialog.titleNode, "inputDialog");
	dojo.addClass(this.dialog.closeButtonNode, "inputDialog");
},

setCombo : function () {

	//while ( this.combo.childNodes )
	//{
	//	this.combo.removeChild(this.combo.childNodes[0]);
	//}
	while ( this.combo.length )
	{
		this.combo.options[this.combo.length - 1] = null;
	}

	for ( var i = 0; i < this.comboValues.length; i++ )
	{
		var option = document.createElement("OPTION");
		option.text = this.comboValues[i];
		option.value = this.comboValues[i];
		this.combo.options.add(option);
	}
},


// LOAD THE DIALOGUE VALUES
load : function (args) {

	this.title 				=	args.title || '';
	this.message 			=	args.message;
	this.comboValues 		=	args.comboValues;
	this.parentWidget 		=	args.parentWidget;
	this.enterCallback 		=	args.enterCallback;
	this.cancelCallback 	=	args.cancelCallback;

	// SET THE DIALOG
	this.dialog.titleNode.innerHTML	=	args.title;
	this.messageNode.innerHTML		=	args.message;
	this.inputMessage.innerHTML		=	args.inputMessage || '';
	this.comboMessage.innerHTML		=	args.comboMessage || '';
	this.checkboxMessageNode.innerHTML	=	args.checkboxMessage || '';
	this.dialog.enterCallback		=	args.enterCallback;
	this.dialog.cancelCallback		=	args.cancelCallback

	// SET ENTER BUTTON AND CANCEL BUTTON LABELS
	this.setEnterLabel(args.enterLabel);
	this.setCancelLabel(args.cancelLabel);

	// SET COMBO BOX IF COMBO MESSAGE IS DEFINED
	if ( this.comboMessage != null )
		this.setCombo();

	this.show();
},

doEnter : function(type) {

	var input = this.inputNode.value;
	var selected = this.combo.value;
	var checked = this.checkbox.checked;
	if ( checked == true ) checked = 1;
	else checked = 0;


	// DO CALLBACK
	this.dialog.enterCallback(input, selected, checked, this);	
}


});
