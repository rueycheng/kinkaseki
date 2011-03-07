#ifndef EXPATPP_CALLBACK_H
#define EXPATPP_CALLBACK_H

#include "expatpp.h"
#include <iostream>

namespace expatpp {
    namespace callback {
	//--------------------------------------------------
	// Callback placeholders
	//-------------------------------------------------- 
	template<class T> void startElementHandler(void* userData, const Char* name, const Char** attr) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->startElement(name, attr);
	}

	template<class T> void endElementHandler(void* userData, const Char* name) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->endElement(name);
	}

	template<class T> void characterDataHandler(void* userData, const Char* s, int len) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->characterData(s, len);
	}

	template<class T> void processingInstructionHandler(void* userData, const Char* target, const Char* data) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->processingInstruction(target, data);
	}

	template<class T> void commentHandler(void* userData, const Char* data) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->comment(userData, data);
	}

	template<class T> void startCdataSectionHandler(void* userData) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->startCdataSection();
	}

	template<class T> void endCdataSectionHandler(void* userData) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->endCdataSection();
	}

	template<class T> void defaultHandler(void* userData, const Char* s, int len) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->handler(s, len);
	}

	template<class T> void skippedEntityHandler(void* userData, const Char* entityName, int isParameterEntity) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->skippedEntity(entityName, isParameterEntity);
	}

	template<class T> void startNamespaceDeclHandler(void* userData, const Char* prefix, const Char* uri) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->startNamespaceDecl(prefix, uri);
	}

	template<class T> void endNamespaceDeclHandler(void* userData, const Char* prefix) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->endNamespaceDecl(prefix);
	}

	template<class T> void xmlDeclHandler(void* userData, const Char* version, const Char* encoding, int standalone) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->xmlDecl(version, encoding, standalone);
	}

	template<class T> void startDoctypeDeclHandler(void* userData, const Char* doctypeName, const Char* sysid,
		const Char* pubid, int hasInternalSubset) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->startDoctypeDecl(doctypeName, sysid, pubid, hasInternalSubset);
	}

	template<class T> void endDoctypeDeclHandler(void* userData) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->endDoctypeDecl();
	}

	template<class T> void elementDeclHandler(void* userData, const Char* name, Content* model) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->elementDecl(name, model);
	}

	template<class T> void attlistDeclHandler(void* userData, const Char* elname, const Char* attname,
		const Char* attType, const Char* dflt, int isRequired) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->attlistDecl(elname, attname, attType, dflt, isRequired);
	}


	template<class T> void entityDeclHandler(void* userData, const Char* entityName, int isParameterEntity,
		const Char* value, int valueLength, const Char* base, const Char* systemId,
		const Char* publicId, const Char* notationName) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->entityDecl(entityName, isParameterEntity, value, valueLength, base, systemId, 
		    publicId, notationName);
	}

	template<class T> void unparsedEntityDeclHandler(void* userData, const Char* entityName, const Char* base,
		const Char* systemId, const Char* publicId, const Char* notationName) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->unparsedEntityDecl(entityName, base, systemId, publicId, notationName);
	}

	template<class T> void notationDeclHandler(void* userData, const Char* notationName, const Char* base,
		const Char* systemId, const Char* publicId) {
	    T* self = static_cast<T*>(userData);
	    if (self) self->notationDecl(notationName, base, systemId, publicId);
	}

	//--------------------------------------------------
	// SAXHandler
	//-------------------------------------------------- 
	class SAXHandler {
	protected:
	    Parser parser;

	public:
	    const static int BUFSIZE = 1 * 1024;
	    SAXHandler() {
		parser.setUserData(static_cast<void*>(this));
		parser.setStartElementHandler(&startElementHandler<SAXHandler>);
		parser.setEndElementHandler(&endElementHandler<SAXHandler>);
		parser.setCharacterDataHandler(&characterDataHandler<SAXHandler>);
	    }

	    virtual void startElement(const Char* name, const Char** attr) {}
	    virtual void endElement(const Char* name) {}
	    virtual void characterData(const Char* s, int len) {}

	    bool carp(std::ostream& err) { 
		err << parser.getCurrentLineNumber() << ',' << parser.getCurrentColumnNumber() << ':' << 
		    errorString(parser.getErrorCode()) << "(" << parser.getCurrentByteCount() << ")\n"; 
		throw 0; // HACK! For debugging purpose
		return true;
	    }

	    bool carp(std::ostream& err, char* buf) { 
		err << parser.getCurrentLineNumber() << ',' << parser.getCurrentColumnNumber() << ':' << 
		    errorString(parser.getErrorCode()) << "(" << parser.getCurrentByteCount() << ")\n"; 
		err << "buffer: \n";
		err << std::string(buf, BUFSIZE) << "\n";
		throw 0; // HACK! For debugging purpose
		return true;
	    }

	    void process(std::istream& in, std::ostream& err = std::cerr) {
		char buf[BUFSIZE];
		while (in.read(buf, BUFSIZE))
		    parser.parse(buf, in.gcount(), false) || carp(err);
		parser.parse(buf, in.gcount(), true) || carp(err);
	    }

	    void process(std::istream& in, const std::string& root = "root", std::ostream& err = std::cerr) {
		char buf[BUFSIZE];
		std::string start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><" + root + ">";
		std::string end = "</" + root + ">";

		parser.parse(start.c_str(), start.size(), false) || carp(err, buf);
		while (in.read(buf, BUFSIZE))
		    parser.parse(buf, in.gcount(), false) || carp(err, buf);
		parser.parse(buf, in.gcount(), false) || carp(err, buf);
		parser.parse(end.c_str(), end.size(), true) || carp(err, buf);
	    }
	};
    }
}
#endif
