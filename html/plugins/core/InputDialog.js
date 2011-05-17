dojo.provide( "plugins.core.InputDialog");

// HAS A
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");


dojo.declare( "plugins.core.InputDialog",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "core/templates/inputdialog.html"),

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins.core") + "/css/inputdialog.css" ],

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT plugins.workflow.Apps WIDGET
parentWidget : null,

// APPLICATION OBJECT
application : null,

// DIALOG TITLE
title: null,

// DISPLAYED MESSAGE 
message : null,

constructor : function(args) {

	this.title 				=	args.title;
	this.message 			=	args.message;
	this.parentWidget 		=	args.parentWidget;
	this.enterCallback 		=	args.enterCallback;
	this.cancelCallback 	=	args.cancelCallback;
	this.checkboxMessage 	=	args.checkboxMessage;

	// SET ENTER BUTTON AND CANCEL BUTTON LABELS
	this.setEnterLabel(args.enterLabel);
	this.setCancelLabel(args.cancelLabel);

	// LOAD CSS
	this.loadCSS();
},

postCreate : function() {

	this.startup();
},

startup : function () {

	this.inherited(arguments);

	// SET UP DIALOG
	this.setDialogue();

	// SHOW CHECKBOX IF this.checkboxMessage IS DEFINED
	this.showCheckbox(this.checkboxMessage);

	// ADD CSS NAMESPACE CLASS FOR TITLE CSS STYLING
	this.setNamespaceClass();
},

// ADD CSS NAMESPACE CLASS
setNamespaceClass : function () {
	dojo.addClass(this.dialog.containerNode, "inputDialog");
	dojo.addClass(this.dialog.titleNode, "inputDialog");
	dojo.addClass(this.dialog.closeButtonNode, "inputDialog");	
},

// SHOW THE DIALOGUE
show: function () {
	this.dialog.show();
},

// HIDE THE DIALOGUE
hide: function () {
	this.dialog.hide();
},

// RUN ENTER CALLBACK IF 'ENTER' CLICKED
doEnter : function(type) {

	var input = this.inputNode.value;
	var checked = this.checkbox.checked;
	if ( checked == true ) checked = 1;
	else checked = 0;


	// DO CALLBACK
	this.dialog.enterCallback(input, checked);	

	// REMOVE INPUT AND HIDE
	this.inputNode.value = '';
	this.dialog.hide();
},

// RUN CANCEL CALLBACK IF 'CANCEL' CLICKED
doCancel : function() {
	this.dialog.cancelCallback();
	this.dialog.hide();
},

// APPEND TO DOCUMENT BODY
setDialogue : function () {

	// APPEND DIALOG TO DOCUMENT
	document.body.appendChild(this.dialog.domNode);
	//this.dialog.show();
},

// LOAD THE DIALOGUE VALUES
load : function (args) {
	//console.dir(this.dialog);

	// SET THE DIALOG
	if ( args.title == null )	{	args.title = "";	}
	this.dialog.titleNode.innerHTML	=	args.title;
	this.messageNode.innerHTML		=	args.message;
	this.dialog.enterCallback		=	args.enterCallback;
	this.dialog.cancelCallback		=	args.cancelCallback

	// SET CHECKBOX BOX IF CHECKBOX MESSAGE IS DEFINED
	this.showCheckbox(args.checkboxMessage);

	// SET ENTER BUTTON AND CANCEL BUTTON LABELS
	this.setEnterLabel(args.enterLabel);
	this.setCancelLabel(args.cancelLabel);

	this.show();
},


// SHOW CHECKBOX
showCheckbox : function (message) {

	if ( message != null && message != '' )
	{
		this.checkboxMessageNode.style.visibility = "visible";
		this.checkbox.style.visibility = "visible";
		this.checkboxMessageNode.innerHTML = message;
	}
	else {
		this.checkboxMessageNode.style.visibility = "hidden";
		this.checkbox.style.visibility = "hidden";
	}
	dojo.removeClass(this.checkboxContainer, "hidden");

},

setEnterLabel : function (label) {

	if ( label != null && this.enterButton != null )
		this.enterButton.set('label', label);	
},

setCancelLabel : function (label) {

	if ( label != null && this.cancelButton != null )
		this.cancelButton.set('label', label);	
}



});

