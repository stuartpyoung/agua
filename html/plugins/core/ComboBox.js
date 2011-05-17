dojo.provide("plugins.core.ComboBox");

// MODIFICATION OF dijit.form.ComboBox TO ALLOW STYLING OF THE
// DROPDOWN COMBO BOX POPUP:
//
// 1. DEFAULT STYLES ARE AS IN dijit.form.ComboBox.
//
// 2. ADJUST STYLES PROGRAMMATICALLY USING:
//
// comboBox.popupClass = "myClass dijitReset dijitMenu"
// 
// ADJUST HEIGHT OF MENU ITEMS TO MATCH YOUR STYLE:
//
// comboBox.itemHeight = 25;
//
//
//	E.G.:
//
//		this.groupCombo.popupClass = "groupusers groupCombo dijitReset dijitMenu";
//		this.groupCombo.wrapperClass = "groupusers dijitPopup";
//		this.groupCombo.itemHeight = 30;
//
// 3. ADD THE FOLLOWING CSS:
// 
//	/* NOTE: THIS IS NECESSARY FOR THE NEXT RULE TO WORK!!! */
// .groupusers .groupCombo .dijitReset .dijitMenu li,
// .groupusers .groupCombo .dijitReset .dijitMenuItem li
// {
// 	/*color: #0FF !important;*/
// }
// .groupusers .dijitMenu li,
// .groupusers .dijitMenuItem li
// {
// 	color: #224 !important;
// 	height: 25px;
// }




dojo.require("dijit.form.ComboBox");

dojo.declare(
	"plugins.core.ComboBoxMixin",
	[ dijit.form.ComboBoxMixin ],
{

	_showResultList: function(){
		this._hideResultList();
		// hide the tooltip
		this.displayMessage("");

		// Position the list and if it's too big to fit on the screen then
		// size it to the maximum possible height
		// Our dear friend IE doesnt take max-height so we need to
		// calculate that on our own every time

		// TODO: want to redo this, see
		//		http://trac.dojotoolkit.org/ticket/3272
		//	and
		//		http://trac.dojotoolkit.org/ticket/4108


		// natural size of the list has changed, so erase old
		// width/height settings, which were hardcoded in a previous
		// call to this function (via dojo.marginBox() call)
		dojo.style(this._popupWidget.domNode, {width: "", height: ""});

		var best = this.open();
		// #3212:
		//		only set auto scroll bars if necessary prevents issues with
		//		scroll bars appearing when they shouldn't when node is made
		//		wider (fractional pixels cause this)
		var popupbox = dojo.marginBox(this._popupWidget.domNode);
		this._popupWidget.domNode.style.overflow =
			((best.h == popupbox.h) && (best.w == popupbox.w)) ? "hidden" : "auto";
		// #4134:
		//		borrow TextArea scrollbar test so content isn't covered by
		//		scrollbar and horizontal scrollbar doesn't appear
		var newwidth = best.w;
		if(best.h < this._popupWidget.domNode.scrollHeight){
			newwidth += 16;
		}
		dojo.marginBox(this._popupWidget.domNode, {
			h: best.h,
			w: Math.max(newwidth, this.domNode.offsetWidth)
		});

		// If we increased the width of drop down to match the width of ComboBox.domNode,
		// then need to reposition the drop down (wrapper) so (all of) the drop down still
		// appears underneath the ComboBox.domNode
		if(newwidth < this.domNode.offsetWidth){
			this._popupWidget.domNode.parentNode.style.left = dojo.position(this.domNode).x + "px";
		}





		// NEED TO CHANGE CLASS OF POPUP WRAPPER TO MAKE POPUP STYLE CHANGES WORK
		var wrapperClass = "dijitPopup";
		if ( this.wrapperClass )	wrapperClass = this.wrapperClass;

		this._popupWidget.domNode.parentNode.setAttribute('class', wrapperClass);







		dijit.setWaiState(this.comboNode, "expanded", "true");
	},


	_onBlur: function(){
		// summary: called magically when focus has shifted away from this widget and it's dropdown


		this._hideResultList();
		//this._arrowIdle();
		//this.inherited(arguments);
	},

});


// INHERITS:
//dojo.declare(
//	"dijit.form._ComboBoxMenu",
//	[dijit._Widget, dijit._Templated],

dojo.declare("plugins.core._ComboBoxMenu",
	[ dijit.form._ComboBoxMenu ],
{
	parentWidget : null,
	popupClass : null,

	constructor : function (args)
	{
		// fill in template with i18n messages

		this.parentWidget = args.parentWidget;
		this.popupClass = args.popupClass;
		this.wrapperClass = args.wrapperClass;
		this.itemHeight = args.itemHeight;
	},


	postCreate:function(){
		// fill in template with i18n messages

		//

		// APPEND TO PARENT WIDGET DOMNODE
		this.parentWidget.domNode.appendChild(this.domNode);

		this.previousButton.innerHTML = this._messages["previousMessage"];
		this.nextButton.innerHTML = this._messages["nextMessage"];
		this.inherited(arguments);
	}

});


dojo.declare("plugins.core.ComboBox",
	[ dijit.form.ComboBox, plugins.core.ComboBoxMixin ],
{

	//postMixInProperties: function(){
	//	// this.inherited(arguments); // ??
	//	dijit.form.ComboBoxMixin.prototype.postMixInProperties.apply(this, arguments);
	//	dijit.form.ValidationTextBox.prototype.postMixInProperties.apply(this, arguments);
	//},

	//
	//_setDisabledAttr: function(/*Boolean*/ value){
	//	dijit.form.ValidationTextBox.prototype._setDisabledAttr.apply(this, arguments);
	//	dijit.form.ComboBoxMixin.prototype._setDisabledAttr.apply(this, arguments);
	//}

});



