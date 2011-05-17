dojo.provide( "plugins.core.CopyWorkflowDialog");

/* CLASS SUMMARY: AN INTERACTIVE DIALOG FOR COPYING WORKFLOWS

	LIKE ITS INHERITED CLASS, InteractiveDialog, CopyWorkflowDialog

	HANGS AROUND UNTIL THE enterCallback METHOD CLOSES IT. IN

	ADDITION, CopyWorkflowDialog ALLOWS THE USER TO SELECT THE

	DESTINATION PROJECT AND THE NAME OF THE NEW WORKFLOW.
*/

// HAS A
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");
dojo.require("plugins.core.InteractiveDialog");


dojo.declare( "plugins.core.CopyWorkflowDialog",
	[ plugins.core.InteractiveDialog ],
{
	//////}}

// SHOW THE DIALOGUE
show: function () {
	this.dialog.show();
},

// HIDE THE DIALOGUE
hide: function () {
	this.dialog.hide();
},

constructor : function(args) {

	this.title 				=	args.title;
	this.message 			=	args.message;
	this.selectValues 		=	args.selectValues;
	this.parentWidget 		=	args.parentWidget;
	this.enterCallback 		=	args.enterCallback;
	this.cancelCallback 	=	args.cancelCallback;

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

	// ADD CSS NAMESPACE CLASS
	dojo.addClass(this.dialog.containerNode, "inputDialog");
	dojo.addClass(this.dialog.titleNode, "inputDialog");
	dojo.addClass(this.dialog.closeButtonNode, "inputDialog");
},


doEnter : function(type) {

	var input = this.inputNode.value;
	// DO CALLBACK
	this.dialog.enterCallback(input, this);

},

// RUN CANCEL CALLBACK IF 'CANCEL' CLICKED
doCancel : function() {
	this.dialog.cancelCallback();
	this.dialog.hide();
},

close : function () {
		// REMOVE INPUT AND HIDE
	this.inputNode.value = '';
	this.dialog.hide();
}



});

