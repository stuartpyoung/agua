dojo.provide( "plugins.core.InteractiveDialog");
/* CLASS SUMMARY: AN INTERACTIVE INPUT DIALOG

	INTERACTIVELY RUN 'Enter' CALLBACK UNTIL CLOSED BY CALLBACK.

	UNLIKE IT'S PARENT CLASS inputDialog, interactiveDialog DOES

	NOT IMMEDIATELY DISAPPEAR AFTER 'Enter' HAS BEEN CLICKED.

	RATHER, IT HANGS AROUND UNTIL THE enterCallback METHOD CLOSES IT.
*/

// HAS A
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");
dojo.require("plugins.core.InputDialog");


dojo.declare( "plugins.core.InteractiveDialog",
	[ plugins.core.InputDialog ],
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

doEnter : function(type) {

	var input = this.inputNode.value;
	var checked = this.checkbox.checked;
	if ( checked == true ) checked = 1;
	else checked = 0;


	// DO CALLBACK
	this.dialog.enterCallback(input, checked, this);		
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

