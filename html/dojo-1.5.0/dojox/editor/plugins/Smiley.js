dojo.provide("dojox.editor.plugins.Smiley");
dojo.experimental("dojox.editor.plugins.Smiley");

dojo.require("dojo.i18n");
dojo.require("dijit._editor._Plugin");
dojo.require("dijit.form.ToggleButton");
dojo.require("dijit.form.DropDownButton");
dojo.require("dojox.editor.plugins._SmileyPalette");

dojo.requireLocalization("dojox.editor.plugins", "Smiley");

dojo.declare("dojox.editor.plugins.Smiley", dijit._editor._Plugin, {
	// summary:
	//		This plugin allows the user to select from emoticons or "smileys"
	//		to insert at the current cursor position.
	//
	// description:
	//		The commands provided by this plugin are:
	//		* smiley - inserts the selected emoticon

	// iconClassPrefix: [const] String
	//		The CSS class name for the button node is formed from `iconClassPrefix` and `command`
	iconClassPrefix: "dijitAdditionalEditorIcon",

	// a marker for emoticon wrap like [:-)] for regexp convienent
	// when a message which contains an emoticon stored in a database or view source, this marker include also
	// but when user enter an emoticon by key board, user don't need to enter this marker.
	// also emoticon definition character set can not contains this marker
	emoticonMarker: '[]',

	emoticonImageClass: 'dojoEditorEmoticon',

	_initButton: function(){
		// summary:
		//
		this.dropDown = new dojox.editor.plugins._SmileyPalette();
		this.connect(this.dropDown, "onChange", function(ascii){
			this.button.closeDropDown();
			this.editor.focus();
			//
			ascii = this.emoticonMarker.charAt(0) + ascii + this.emoticonMarker.charAt(1);
			this.editor.execCommand("inserthtml", ascii);
		});
		this.i18n = dojo.i18n.getLocalization("dojox.editor.plugins", "Smiley");
		this.button = new dijit.form.DropDownButton({
			label: this.i18n.smiley,
			showLabel: false,
			iconClass: this.iconClassPrefix + " " + this.iconClassPrefix + "Smiley",
			tabIndex: "-1",
			dropDown: this.dropDown
		});
		this.emoticonImageRegexp = new RegExp("class=(\"|\')" + this.emoticonImageClass + "(\"|\')");
	},

	setEditor: function(editor){
		// summary:
		//		Over-ride for the setting of the editor.
		// editor: Object
		//		The editor to configure for this plugin to use.
		this.editor = editor;
		this._initButton();
		this.editor.contentPreFilters.push(dojo.hitch(this, this._preFilterEntities));
		this.editor.contentPostFilters.push(dojo.hitch(this, this._postFilterEntities));
	},

	_preFilterEntities: function(/*String content passed in*/ value){
		// summary:
		//		A function to filter out emoticons into their UTF-8 character form
		//		displayed in the editor.  It gets registered with the preFilters
		//		of the editor.
		// tags:
		//		private.
		//
		//
		return value.replace(/\[([^\]]*)\]/g, dojo.hitch(this, this._decode));
	},

	_postFilterEntities: function(/*String content passed in*/ value){
		// summary:
		//		A function to filter out emoticons into encoded form so they
		//		are properly displayed in the editor.  It gets registered with the
		//		postFilters of the editor.
		// tags:
		//		private.
		return value.replace(/<img [^>]*>/gi, dojo.hitch(this, this._encode));
	},

	_decode: function(str, ascii){
		// summary:
		//		Pre-filter for editor to convert strings like [:-)] into an <img> of the corresponding smiley
		var emoticon = dojox.editor.plugins.Emoticon.fromAscii(ascii);
		return emoticon ? emoticon.imgHtml(this.emoticonImageClass) : ascii;
	},

	_encode: function(str){
		// summary:
		//		Post-filter for editor to convert <img> nodes of smileys into strings like [:-)]

		// Each <img> node has an alt tag with it's ascii representation, so just use that.
		// TODO: wouldn't this be easier as a postDomFilter ?
		if(str.search(this.emoticonImageRegexp) > -1){
			return this.emoticonMarker.charAt(0) + str.replace(/(<img [^>]*)alt="([^"]*)"([^>]*>)/, "$2") + this.emoticonMarker.charAt(1);
		}
		else{
			return str;
		}
	}
});

// Register this plugin.
dojo.subscribe(dijit._scopeName + ".Editor.getPlugin",null,function(o){
	if(o.plugin){ return; }
	if(o.args.name === "smiley"){
		o.plugin = new dojox.editor.plugins.Smiley();
	}
});
