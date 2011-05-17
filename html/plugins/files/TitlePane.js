
dojo.provide("plugins.files.TitlePane");

// INHERITS
dojo.require("dijit.TitlePane");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.files.FileDrag");

dojo.declare(
	"plugins.files.TitlePane",
	[dijit.TitlePane, plugins.core.Common],
{
	// summary: A pane with a title on top, that can be opened or collapsed.

	// className : String. Name of class
	className: "filesTitlePane",

	// title: String. Title of the pane
	title: "",

	// open: Boolean. Whether pane is opened or closed.
	open: true,

	// duration: Integer
	//		Time in milliseconds to fade in/fade out
	duration: dijit.defaultDuration,

	// baseClass: String
	//	The root className to use for the various states of this widget
	baseClass: "filesTitlePane",

	// templateString: String, set to null so that template is loaded from templatePath
	templateString : null,

	// templatePath: String, Path to html template
	templateString:"<div dojoAttachPoint=\"containerNode\">\r\n\r\n\t\t<div dojoAttachEvent=\"onclick:toggle,onkeypress: _onTitleKey,onfocus:_handleFocus,onblur:_handleFocus\" tabindex=\"0\"\r\n\t\t\twaiRole=\"button\"\r\n\t\t\tclass=\"dijitTitlePaneTitle\" dojoAttachPoint=\"titleBarNode,focusNode\">\r\n\t\r\n\t\t\t<img src=\"${_blankGif}\" alt=\"\" dojoAttachPoint=\"arrowNode\" class=\"filesArrowNode\" waiRole=\"presentation\">\r\n\t\t\t\t<span dojoAttachPoint=\"arrowNodeInner\" class=\"dijitArrowNodeInner\"></span>\r\n\t\t\t\t<img src=\"${_blankGif}\" dojoAttachPoint=\"refresh\" class=\"refreshTitlePane\" dojoAttachEvent=\"onclick:reload\">\r\n\t\r\n\t\t\t<span dojoAttachPoint=\"titleNode\" class=\"dijitTitlePaneTextNode\">\r\n\t\t\t\t\r\n\t\t\t\t<table dojoAttachPoint='header' class='header'>\r\n\t\t\t\t\t<tr>\r\n\t\t\t\t\t\t<td dojoAttachPoint=\"ownerNode\" class='owner' title='Owner'></td>\r\n\t\t\t\t\t\t<td dojoAttachPoint=\"nameNode\" class='name'></td>\r\n\t\t\t\t\t\t<td dojoAttachPoint=\"descriptionNode\" class='description' title='Description'></td>\r\n\t\t\t\t\t</tr>\r\n\t\t\t\t</table>\r\n\t\t\t\t\r\n\t\t\t</span>\r\n\t\t\t\r\n\t\t</div>\r\n\t\t<div class=\"dijitTitlePaneContentOuter\" dojoAttachPoint=\"hideNode\">\r\n\t\t\t<div class=\"dijitReset\" dojoAttachPoint=\"wipeNode\">\r\n\t\t\t\t<div class=\"dijitTitlePaneContentInner\" dojoAttachPoint=\"containerNode\" waiRole=\"region\" tabindex=\"-1\">\r\n\t\t\t\t\t<!-- nested divs because wipeIn()/wipeOut() doesn't work right on node w/padding etc.  Put padding on inner div. -->\r\n\t\t\t\t</div>\r\n\t\t\t</div>\r\n\t\t</div>\r\n</div>\r\n",

	// ccsFiles: Array. Paths to CSS files
	cssFiles: [ "plugins/files/css/titlepane.css" ],

	// reloadCallback : Function. Function to call when 'reload' button clicked
	reloadCallback : null,

	// parentWidget: Object. The XXXFiles.js object containing this file title pane
	parentWidget : null,

	// directory: String. Filesystem location
	directory: "",

	// name: String. Name of filesystem
	name: "",

	//// ownerNode: DIV element, insert owner name into its innerHTML
	//ownerNode : null,

	// size: Integer. Size of icons and fonts in title
	size: 'normal',

	// owner: String. Owner of filesystem
	owner: "",

	// type: String. Type of filesystem
	type: "",

	// description: String. Description of filesystem
	description: "",

	attributeMap: dojo.mixin(dojo.clone(dijit.layout.ContentPane.prototype.attributeMap), {
		title: {node: "titleNode", type: "innerHTML" }
	}),

	// RELOAD FileDrag OBJECT INSIDE TITLE PANE	
	reload : function (event)
	{
		event.stopPropagation();


		if ( this.reloadCallback != null )
		{
			this.reloadCallback();
			return;
		}
		//console.dir(this);
		//if ( this.reloadCallback == null )	return;

		// REMOVE EXISTING FILE DRAG
		while ( this.containerNode.firstChild )
		{
			this.containerNode.removeChild(this.containerNode.firstChild);
		}

		var fileDrag = this.parentWidget.createFileDrag(this.directory);

		// ADD fileDrag TO TITLE PANE
		this.containerNode.appendChild(fileDrag.domNode);

		//var childNode = this.containerNode.firstChild;
		//var oldFileDrag = dijit.byNode(childNode);	

		//this.reloadCallback('reloading');
	},


	constructor: function(args){
		//console.dir(this);
		//

		this.loadCSS();

	},

	//
	//
	////Inherited from dijit._Widget and called just before template
	////instantiation in buildRendering. This method is especially useful
	////for manipulating the template before it becomes visible.
	//postMixInProperties: function(args)
	//{
	//	//console.dir(this);
	//
	//},
	//
	//
	postCreate: function(){


		this.loadCSS(this.cssFiles);



		if ( this.owner != null ) this.ownerNode.innerHTML = this.owner;
		if ( this.name != null ) this.nameNode.innerHTML = this.name;
		if ( this.description != null ) this.descriptionNode.innerHTML = this.description;

		// SET TITLE FOR FILESYSTEM TYPE
		var nameTitle = "Type of filesystem";
		//if ( this.type != null ) nameTitle = this.type;
		if ( this.name != null ) this.nameNode.setAttribute('title', nameTitle);

		// ADJUST SIZE IF SPECIFIED

		if ( this.size == 'large' )
		{
			this.arrowNode.setAttribute('class', 'largeFilesArrowNode');

		}

		if ( this.size == 'larger' )
		{
			this.arrowNode.setAttribute('class', 'largerFilesArrowNode');
		}

		this.inherited(arguments);
	}

	//,
	//
	//_setOpenAttr: function(/* Boolean */ open){
	//	// summary:
	//	//		Hook to make attr("open", boolean) control the open/closed state of the pane.
	//	// open: Boolean
	//	//		True if you want to open the pane, false if you want to close it.
	//	if(this.open !== open){ this.toggle(); }
	//},
	//
	//_setContentAttr: function(content){
	//	// summary:
	//	//		Hook to make attr("content", ...) work.
	//	// 		Typically called when an href is loaded.  Our job is to make the animation smooth
	//
	//	if(!this.open || !this._wipeOut || this._wipeOut.status() == "playing"){
	//		// we are currently *closing* the pane (or the pane is closed), so just let that continue
	//		this.inherited(arguments);
	//	}else{
	//		if(this._wipeIn && this._wipeIn.status() == "playing"){
	//			this._wipeIn.stop();
	//		}
	//
	//		// freeze container at current height so that adding new content doesn't make it jump
	//		dojo.marginBox(this.wipeNode, { h: dojo.marginBox(this.wipeNode).h });
	//
	//		// add the new content (erasing the old content, if any)
	//		this.inherited(arguments);
	//
	//		// call _wipeIn.play() to animate from current height to new height
	//		if(this._wipeIn){
	//			this._wipeIn.play();
	//		}else{
	//			this.hideNode.style.display = "";
	//		}
	//	}
	//},

//	toggle : function(){
//
//
//		// summary: switches between opened and closed state
//		dojo.forEach([this._wipeIn, this._wipeOut], function(animation){
//			if(animation && animation.status() == "playing"){
//				animation.stop();
//			}
//		});
//
//		var anim = this[this.open ? "_wipeOut" : "_wipeIn"]
//		if(anim){
//			anim.play();
//		}else{
//			this._hideNode.style.display = this.open ? "" : "none";
//		}
//		this.open =! this.open;
//
//		// load content (if this is the first time we are opening the TitlePane
//		// and content is specified as an href, or href was set when hidden)
//		this._loadCheck();
//
//		this._setCss();
//	},

	//_setCss: function(){
	//	// summary: set the open/close css state for the TitlePane
	//	var classes = ["dijitClosed", "dijitOpen"];
	//	var boolIndex = this.open;
	//	var node = this.titleBarNode || this.focusNode;
	//	dojo.removeClass(node, classes[!boolIndex+0]);
	//	node.className += " " + classes[boolIndex+0];
	//
	//	// provide a character based indicator for images-off mode
	//	this.arrowNodeInner.innerHTML = this.open ? "-" : "+";
	//},
	//
	//_onTitleKey: function(/*Event*/ e){
	//	// summary: callback when user hits a key
	//	if(e.charOrCode == dojo.keys.ENTER || e.charOrCode == ' '){
	//		this.toggle();
	//	}else if(e.charOrCode == dojo.keys.DOWN_ARROW && this.open){
	//		this.containerNode.focus();
	//		e.preventDefault();
	// 	}
	//},
	//

	,
	_handleFocus: function(/*Event*/ e){
		// summary: handle blur and focus for this widget

		// add/removeClass is safe to call without hasClass in this case
		dojo[(e.type == "focus" ? "addClass" : "removeClass")](this.focusNode, this.baseClass + "Focused");
	},
	//
	//setTitle: function(/*String*/ title){
	//	// summary: sets the text of the title
	//	dojo.deprecated("TitlePane.setTitle() is deprecated.  Use attr('title', ...) instead.", "", "2.0");
	//	this.titleNode.innerHTML = title;
	//}

});
