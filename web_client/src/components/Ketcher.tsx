import { useEffect, useRef, useState, useImperativeHandle, forwardRef } from 'react';
import { Box, LoadingOverlay, Text } from '@mantine/core';

export interface KetcherRef {
    getSmiles: () => Promise<string>;
    setMolecule: (molfile: string) => void;
}

interface KetcherProps {
    onStructChange?: (smiles: string) => void;
    initialMolfile?: string;
}

export const Ketcher = forwardRef<KetcherRef, KetcherProps>(({ onStructChange, initialMolfile }, ref) => {
    const iframeRef = useRef<HTMLIFrameElement>(null);
    const [isReady, setIsReady] = useState(false);
    const [ketcherInstance, setKetcherInstance] = useState<any>(null);

    useImperativeHandle(ref, () => ({
        getSmiles: async () => {
            if (!ketcherInstance) return '';
            return await ketcherInstance.getSmiles();
        },
        setMolecule: (molfile: string) => {
            if (ketcherInstance) {
                ketcherInstance.setMolecule(molfile);
            }
        }
    }));

    // Initialize Ketcher connection
    useEffect(() => {
        const iframe = iframeRef.current;
        if (!iframe) return;

        const checkKetcher = () => {
            if (iframe.contentWindow && (iframe.contentWindow as any).ketcher) {
                const ketcher = (iframe.contentWindow as any).ketcher;
                setKetcherInstance(ketcher);
                setIsReady(true);
            } else {
                setTimeout(checkKetcher, 500);
            }
        };

        iframe.onload = checkKetcher;
        checkKetcher(); 
    }, []);

    // Handle initial molfile
    useEffect(() => {
        if (isReady && ketcherInstance && initialMolfile) {
            ketcherInstance.setMolecule(initialMolfile);
        }
    }, [isReady, ketcherInstance, initialMolfile]);
    
    // Fragment detection watcher
    useEffect(() => {
        if (ketcherInstance && onStructChange) {
             let lastSmiles = '';
             const interval = setInterval(async () => {
                 try {
                     const smiles = await ketcherInstance.getSmiles();
                     if (smiles !== lastSmiles) {
                         lastSmiles = smiles;
                         onStructChange(smiles);
                     }
                 } catch (e) {
                     // ignore error
                 }
             }, 1500); // Check every 1.5s
             return () => clearInterval(interval);
        }
    }, [ketcherInstance, onStructChange]);

    return (
        <Box pos="relative" h={500} w="100%" style={{ border: '1px solid #dee2e6', borderRadius: '4px' }}>
            <LoadingOverlay visible={!isReady} />
            {!isReady && <Text ta="center" mt="md">Loading Ketcher...</Text>}
            <iframe
                ref={iframeRef}
                src="/ketcher/index.html?hiddenControls=open"
                width="100%"
                height="100%"
                style={{ border: 'none' }}
                title="Ketcher"
                onLoad={() => window.scrollTo({ top: 0, behavior: 'instant' })}
            />
        </Box>
    );
});
