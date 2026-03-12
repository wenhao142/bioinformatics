import "./globals.css";
import { ReactNode } from "react";

export const metadata = {
  title: "AD Locus Evidence Platform",
  description: "AD multi-omics workflow, locus review, and evidence console",
};

export default function RootLayout({ children }: { children: ReactNode }) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
