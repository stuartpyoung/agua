
dojo.provide("plugins.upload.FileInput");

// FILE UPLOAD AND AUTO UPLOAD
dojo.require("dojox.form.FileInput");	
dojo.require("dojox.form.FileInputAuto");

dojo.declare("plugins.upload.FileInput",
	dojox.form.FileInputAuto,
	{
	// summary: An extended version of dojox.form.FileInputAuto
	//	that does not display an input node, but rather only a button
	// 	and otherwise behaves just like FileInputAuto

	templateString:"<div class=\"dijitFileInput\">\r\n\t<input id=\"${id}\" name=\"${name}\" class=\"dijitFileInputReal\" type=\"file\" dojoAttachPoint=\"fileInput\" />\r\n\t<div class=\"dijitFakeInput\" dojoAttachPoint=\"fakeNodeHolder\">\r\n\t\t<input class=\"dijitFileInputVisible\" type=\"text\" dojoAttachPoint=\"focusNode, inputNode\" />\r\n\t\t<div class=\"dijitInline dijitFileInputText\" dojoAttachPoint=\"titleNode\">${label}</div>\r\n\t\t<div class=\"dijitInline dijitFileInputButton\" dojoAttachPoint=\"cancelNode\" dojoAttachEvent=\"onclick:reset\">${cancelText}</div>\r\n\t</div>\r\n\t<div class=\"dijitProgressOverlay\" dojoAttachPoint=\"overlay\">&nbsp;</div>\r\n</div>\r\n",

	startup: function(){


		// summary: hide our fileInput input field
		this.inherited(arguments);
		this._off = dojo.style(this.inputNode,"width");
		this.inputNode.style.display = "none";
		this._fixPosition();

		this.titleNode
	},

	_handleSend: function(data,ioArgs){
		// summary: The callback to toggle the progressbar, and fire the user-defined callback


		// innerHTML throws errors in IE! so use DOM manipulation instead
		//this.overlay.removeChild(this.overlay.firstChild);

		this._sent = true;
		this._sending = false;

		// DO NOTHING, EXCEPT

		//remove the form used to send the request
		dojo.body().removeChild(ioArgs.args.form);

		this.onComplete(data,ioArgs,this);
	},

	// EVENT WHICH TRIGGERS THE FILE UPLOAD
	triggerEvent : "onchange",

	// KILL onBlur EVENT
	_onBlur: function(){

		// RETURN IF NOT JUST OPENED BROWSER
		if ( this._openedBrowser != true )   return;

		// summary: start the upload timer
		if ( this._blurTimer ) { clearTimeout(this._blurTimer); }
		if ( !this._sent ){
			this._blurTimer = setTimeout(dojo.hitch(this, "_sendFile"),this.blurDelay);		
		}
	},

	reset: function(/* Event */e){


		// DO NOTHING, EXCEPT
		this._fixPosition(); 
	},

	_fixPosition: function(){		


		// summary: in this case, set the button under where the visible button is 
		if(dojo.isIE){
			dojo.style(this.fileInput,"width","1px");
		}else{
			dojo.style(this.fileInput,"left","-"+(this._off)+"px");
		}
	},

	// TRIGGER THE CHAIN OF EVENTS TO UPLOAD A FILE
	_sendFile : function(e){
		// IF THE INPUT ELEMENT HAS NO VALUE, ABORT THE SEND AND DESTROY THIS OBJECT

		if (this._sent || this._sending || !this.fileInput.value || this._openedBrowser == false )
		{
			return;
		}
		this._sending = true;
		this._openedBrowser = false; 

		this.overlay.setAttribute('style', 'visibility: visible; height: 20px; width: auto;');

		dojo.style(this.fakeNodeHolder,"display","none");
		dojo.style(this.overlay,{
			opacity:0,
			display:"block"
		});

		this.setMessage(this.uploadMessage);
		dojo.fadeIn({ node: this.overlay, duration: this.duration }).play();

		var _newForm; 
		if(dojo.isIE){
			// just to reiterate, IE is a steaming pile of code. 
			_newForm = document.createElement('<form enctype="multipart/form-data" method="post">');
			_newForm.encoding = "multipart/form-data";	
		}
		else
		{
			// this is how all other sane browsers do it
			_newForm = document.createElement('form');
			_newForm.setAttribute("enctype","multipart/form-data");
		}

		// ADD CLONED fileInput TO NEW FORM IN ORDER TO UPLOAD FILE.
		// DON'T ADD THE ORIGINAL OR IT'LL BE DISPLACED AND/OR DESTROYED
		var clonedFileInput = dojo.clone(this.fileInput);
		_newForm.appendChild(clonedFileInput);	
		dojo.body().appendChild(_newForm);

		// AUTOMATED CREATION OF HIDDEN NODES
		var hiddenValues = this.hiddenValues;
		for ( var key in hiddenValues )
		{
			var node = document.createElement('input');
			node.type = "hidden";
			node.name = key;
			node.value = hiddenValues[key];
			node.id = "hidden" + key;
			_newForm.appendChild(node);
		}
		dojo.io.iframe.send({
			url: this.url,
			form: _newForm,
			timeout : 1000,
			handleAs: "json",
			handle: dojo.hitch(this,"onComplete")
		});

		// GET RID OF FORM TO AVOID IT FIRING ON ANY OTHER RANDOM ONBLUR EVENTS
		_newForm.removeChild(clonedFileInput);	
		dojo.body().appendChild(_newForm);
	},


	onComplete : function(){

		// REMOVE 'Uploading file' MESSAGE

		this.overlay.setAttribute('style', 'visibility: hidden;');

		// SET THE UPLOADED FILE NAME IN THE TABLE ELEMENT 

		// SET THIS.SENDING TO FALSE
		this._sending = false;

		// DO CALLBACK TO UPDATE THE WORKFLOW AND THE OPTIONS TABLE
		if ( this.fileInput.value == null || this.fileInput.value == '' )
		{
			this.callback();
		}
		else
		{
			var filepath = this.path + "/" + this.fileInput.value;
			this.callback(filepath);				
		}

	},


	reset : function(/* Event */e){
		// summary: on click of cancel button, since we can't clear the input because of
		// 	security reasons, we destroy it, and add a new one in it's place.



		this.disconnect(this._listener);
		this.disconnect(this._keyListener);
		if(this.fileInput){

			this.domNode.removeChild(this.fileInput);
		}
		dojo.fadeOut({ node: this.cancelNode, duration:275 }).play(); 

		// should we use cloneNode()? can we?
		this.fileInput = document.createElement('input');
		// dojo.attr(this.fileInput,{
		//	"type":"file", "id":this.id, "name": this.name	
		//});
		this.fileInput.setAttribute("FileUpload.generateFileInput    type","file");
		this.fileInput.setAttribute("FileUpload.generateFileInput    id", this.id);
		this.fileInput.setAttribute("FileUpload.generateFileInput    name", this.name);
		dojo.addClass(this.fileInput,"dijitFileInputReal");
		this.domNode.appendChild(this.fileInput);

		this._keyListener = this.connect(this.fileInput, "onkeyup", "_matchValue");
		this._listener = this.connect(this.fileInput, "onchange", "_matchValue"); 
		this.inputNode.value = ""; 
	}



});
